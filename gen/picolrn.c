/*
 *  This code is released under the MIT License.
 *  Copyright (c) 2013 Nenad Markus
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include <omp.h>

// hyperparameters
#define NRANDS 128

/*
	auxiliary stuff
*/

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define SQR(x) ((x)*(x))

/*
	portable time function
*/

#ifdef __GNUC__
#include <time.h>
float getticks()
{
	struct timespec ts;

	if(clock_gettime(CLOCK_MONOTONIC, &ts) < 0)
		return -1.0f;

	return ts.tv_sec + 1e-9f*ts.tv_nsec;
}
#else
#include <windows.h>
float getticks()
{
	static double freq = -1.0;
	LARGE_INTEGER lint;

	if(freq < 0.0)
	{
		if(!QueryPerformanceFrequency(&lint))
			return -1.0f;

		freq = lint.QuadPart;
	}

	if(!QueryPerformanceCounter(&lint))
		return -1.0f;

	return (float)( lint.QuadPart/freq );
}
#endif

/*
	multiply with carry PRNG
*/

uint32_t mwcrand_r(uint64_t* state)
{
	uint32_t* m;

	//
	m = (uint32_t*)state;

	// bad state?
	if(m[0] == 0)
		m[0] = 0xAAAA;

	if(m[1] == 0)
		m[1] = 0xBBBB;

	// mutate state
	m[0] = 36969 * (m[0] & 65535) + (m[0] >> 16);
	m[1] = 18000 * (m[1] & 65535) + (m[1] >> 16);

	// output
	return (m[0] << 16) + m[1];
}

uint64_t prngglobal = 0x12345678000fffffLL;

void smwcrand(uint32_t seed)
{
	prngglobal = 0x12345678000fffffLL*seed;
}

uint32_t mwcrand()
{
	return mwcrand_r(&prngglobal);
}

/*
	
*/

#define MAX_N 4000000

int nimages = 0;
uint8_t* ppixels[MAX_N];
int pdims[MAX_N][2]; // (nrows, ncols)
int contents[MAX_N][2]; // (i, j) where i is the start index and j-1 is the ending index

int nobjects = 0;
int objects[MAX_N][4]; // (r, c, s, i)

int load_image(uint8_t* pixels[], int* nrows, int* ncols, FILE* file)
{
	/*
	- loads an 8-bit grey image saved in the <RID> file format
	- <RID> file contents:
		- a 32-bit signed integer h (image height)
		- a 32-bit signed integer w (image width)
		- an array of w*h unsigned bytes representing pixel intensities
	*/

	//
	if(fread(nrows, sizeof(int), 1, file) != 1)
		return 0;

	if(fread(ncols, sizeof(int), 1, file) != 1)
		return 0;

	//
	*pixels = (uint8_t*)malloc(*nrows**ncols*sizeof(uint8_t));

	if(!*pixels)
		return 0;

	// read pixels
	if(fread(*pixels, sizeof(uint8_t), *nrows**ncols, file) != *nrows**ncols)
		return 0;

	// we're done
	return 1;
}

int load_training_data(char* path)
{
	FILE* file;

	//
	file = fopen(path, "rb");

	if(!file)
		return 0;

	//
	nimages = 0;
	nobjects = 0;

	while( load_image(&ppixels[nimages], &pdims[nimages][0], &pdims[nimages][1], file) )
	{
		int i, n;

		//
		contents[nimages][0] = nobjects;

		if(fread(&n, sizeof(int), 1, file) != 1)
			return 1;

		for(i=0; i<n; ++i)
		{
			fread(&objects[nobjects][0], sizeof(int), 1, file); // r
			fread(&objects[nobjects][1], sizeof(int), 1, file); // c
			fread(&objects[nobjects][2], sizeof(int), 1, file); // s

			objects[nobjects][3] = nimages; // i

			//
			++nobjects;
		}

		contents[nimages][1] = contents[nimages][0] + n;

		//
		++nimages;
	}

	//
	return 1;
}

void delete_training_data()
{
	int i;

	for(i=0; i<nimages; ++i)
	{
		free(ppixels[i]);
		ppixels[i] = 0;
	}

	nimages = 0;
}

/*
	regression trees
*/

int bintest(int32_t tcode, int r, int c, int s, int iind)
{
	//
	int r1, c1, r2, c2;
	int8_t* p = (int8_t*)&tcode;

	//
	r1 = (256*r + p[0]*s)/256;
	c1 = (256*c + p[1]*s)/256;

	r2 = (256*r + p[2]*s)/256;
	c2 = (256*c + p[3]*s)/256;

	//
	r1 = MIN(MAX(0, r1), pdims[iind][0]-1);
	c1 = MIN(MAX(0, c1), pdims[iind][1]-1);

	r2 = MIN(MAX(0, r2), pdims[iind][0]-1);
	c2 = MIN(MAX(0, c2), pdims[iind][1]-1);

	//
	return ppixels[iind][r1*pdims[iind][1]+c1]<=ppixels[iind][r2*pdims[iind][1]+c2];
}

float get_split_error(int32_t tcode, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int inds[], int indsnum)
{
	int i;

	double wsum, wsum0, wsum1;
	double wtvalsum0, wtvalsumsqr0, wtvalsum1, wtvalsumsqr1;

	double wmse0, wmse1;

	//
	wsum = wsum0 = wsum1 = wtvalsum0 = wtvalsum1 = wtvalsumsqr0 = wtvalsumsqr1 = 0.0;

	for(i=0; i<indsnum; ++i)
	{
		if( bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], iinds[inds[i]]) )
		{
			wsum1 += ws[inds[i]];
			wtvalsum1 += ws[inds[i]]*tvals[inds[i]];
			wtvalsumsqr1 += ws[inds[i]]*SQR(tvals[inds[i]]);
		}
		else
		{
			wsum0 += ws[inds[i]];
			wtvalsum0 += ws[inds[i]]*tvals[inds[i]];
			wtvalsumsqr0 += ws[inds[i]]*SQR(tvals[inds[i]]);
		}

		wsum += ws[inds[i]];
	}

	//
	wmse0 = wtvalsumsqr0 - SQR(wtvalsum0)/wsum0;
	wmse1 = wtvalsumsqr1 - SQR(wtvalsum1)/wsum1;

	//
	return (float)( (wmse0 + wmse1)/wsum );
}

int split_training_data(int32_t tcode, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int inds[], int ninds)
{
	int stop;
	int i, j;

	int n0;

	//
	stop = 0;

	i = 0;
	j = ninds - 1;

	while(!stop)
	{
		//
		while( !bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], iinds[inds[i]]) )
		{
			if( i==j )
				break;
			else
				++i;
		}

		while( bintest(tcode, rs[inds[j]], cs[inds[j]], ss[inds[j]], iinds[inds[j]]) )
		{
			if( i==j )
				break;
			else
				--j;
		}

		//
		if( i==j )
			stop = 1;
		else
		{
			// swap
			inds[i] = inds[i] ^ inds[j];
			inds[j] = inds[i] ^ inds[j];
			inds[i] = inds[i] ^ inds[j];
		}
	}

	//
	n0 = 0;

	for(i=0; i<ninds; ++i)
		if( !bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], iinds[inds[i]]) )
			++n0;

	//
	return n0;
}

int32_t get_random_tcode(int8_t* bbox)
{
	int32_t tcode;
	int8_t* p;

	//
	p = (int8_t*)&tcode;

	//
	p[0] = bbox[0] + mwcrand()%(bbox[1]-bbox[0]+1);
	p[1] = bbox[2] + mwcrand()%(bbox[3]-bbox[2]+1);

	p[2] = bbox[0] + mwcrand()%(bbox[1]-bbox[0]+1);
	p[3] = bbox[2] + mwcrand()%(bbox[3]-bbox[2]+1);

	//
	return tcode;
}

int grow_subtree(int32_t tcodes[], float lut[], int nodeidx, int d, int maxd, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int inds[], int ninds, int8_t* bbox)
{
	int i, nrands;

	int32_t tmptcodes[2048];
	float es[2048], e;

	int n0;

	//
	if(d == maxd)
	{
		int lutidx;
		double tvalaccum, wsum;

		//
		lutidx = nodeidx - ((1<<maxd)-1);

		// compute output: a simple average
		tvalaccum = 0.0;
		wsum = 0.0;

		for(i=0; i<ninds; ++i)
		{
			tvalaccum += ws[inds[i]]*tvals[inds[i]];
			wsum += ws[inds[i]];
		}

		if(wsum == 0.0)
			lut[lutidx] = 0.0f;
		else
			lut[lutidx] = (float)( tvalaccum/wsum );

		//
		return 1;
	}
	else if(ninds <= 1)
	{
		//
		tcodes[nodeidx] = 0;

		//
		grow_subtree(tcodes, lut, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, ss, iinds, ws, inds, ninds, bbox);
		grow_subtree(tcodes, lut, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, ss, iinds, ws, inds, ninds, bbox);

		return 1;
	}

	// generate binary test codes
	nrands = NRANDS;

	for(i=0; i<nrands; ++i)
		tmptcodes[i] = get_random_tcode(bbox);

	//
	#pragma omp parallel for
	for(i=0; i<nrands; ++i)
		es[i] = get_split_error(tmptcodes[i], tvals, rs, cs, ss, iinds, ws, inds, ninds);

	//
	e = es[0];
	tcodes[nodeidx] = tmptcodes[0];

	for(i=1; i<nrands; ++i)
		if(e > es[i])
		{
			e = es[i];
			tcodes[nodeidx] = tmptcodes[i];
		}

	//
	n0 = split_training_data(tcodes[nodeidx], tvals, rs, cs, ss, iinds, ws, inds, ninds);

	//
	grow_subtree(tcodes, lut, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, ss, iinds, ws, &inds[0], n0, bbox);
	grow_subtree(tcodes, lut, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, ss, iinds, ws, &inds[n0], ninds-n0, bbox);

	//
	return 1;
}

int grow_rtree(int32_t tcodes[], float lut[], int d, float tvals[], int rs[], int cs[], int ss[], int iinds[], double ws[], int n, int8_t* bbox)
{
	int i;
	int* inds;

	//
	inds = (int*)malloc(n*sizeof(int));

	for(i=0; i<n; ++i)
		inds[i] = i;

	//
	if(!grow_subtree(tcodes, lut, 0, 0, d, tvals, rs, cs, ss, iinds, ws, inds, n, bbox))
	{
		free(inds);
		return 0;
	}
	else
	{
		free(inds);
		return 1;
	}
}

/*
	
*/

int32_t version = 3;

int tdepth;
int ntrees=0;

int8_t bbox[4]; // (r_min, r_max, c_min, c_max)

int32_t tcodes[4096][1024];
float luts[4096][1024];

float thresholds[4096];

/*
	
*/

int load_cascade_from_file(FILE* file)
{
	int i;

	//
	fread(&version, sizeof(int32_t), 1, file);
	fread(&bbox[0], sizeof(int8_t), 4, file);
	fread(&tdepth, sizeof(int), 1, file);
	fread(&ntrees, sizeof(int), 1, file);

	for(i=0; i<ntrees; ++i)
	{
		//
		fread(&tcodes[i][0], sizeof(int32_t), (1<<tdepth)-1, file);
		fread(&luts[i][0], sizeof(float), 1<<tdepth, file);
		fread(&thresholds[i], sizeof(float), 1, file);
	}

	//
	return 1;
}

int save_cascade_to_file(FILE* file)
{
	int i;

	//
	fwrite(&version, sizeof(int32_t), 1, file);
	fwrite(&bbox[0], sizeof(int8_t), 4, file);
	fwrite(&tdepth, sizeof(int), 1, file);
	fwrite(&ntrees, sizeof(int), 1, file);

	for(i=0; i<ntrees; ++i)
	{
		//
		fwrite(&tcodes[i][0], sizeof(int32_t), (1<<tdepth)-1, file);
		fwrite(&luts[i][0], sizeof(float), 1<<tdepth, file);
		fwrite(&thresholds[i], sizeof(float), 1, file);
	}

	//
	return 1;
}

/*
	
*/

float get_tree_output(int i, int r, int c, int s, int iind)
{
	int idx, j;

	//
	idx = 1;

	for(j=0; j<tdepth; ++j)
		idx = 2*idx + bintest(tcodes[i][idx-1], r, c, s, iind);

	//
	return luts[i][idx - (1<<tdepth)];
}

int classify_region(float* o, int r, int c, int s, int iind)
{
	int i;

	//
	*o = 0.0f;

	if(!ntrees)
		return 1;

	//
	for(i=0; i<ntrees; ++i)
	{
		//
		*o += get_tree_output(i, r, c, s, iind);

		//
		if(*o <= thresholds[i])
			return -1;
	}

	//
	return 1;
}

float learn_new_stage(float mintpr, float maxfpr, int maxntrees, float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int np, int nn)
{
	int i;

	double* ws;
	double wsum;

	float threshold, tpr, fpr;

	//
	printf("* learning a new stage ...\n");

	//
	ws = (double*)malloc((np+nn)*sizeof(double));

	//
	maxntrees = ntrees + maxntrees;
	fpr = 1.0f;

	while(ntrees<maxntrees && fpr>maxfpr)
	{
		float t;
		int numtps, numfps;

		//
		t = getticks();

		// compute weights ...
		wsum = 0.0;

		for(i=0; i<np+nn; ++i)
		{
			if(tvals[i] > 0)
				ws[i] = exp(-1.0*os[i])/np;
			else
				ws[i] = exp(+1.0*os[i])/nn;

			wsum += ws[i];
		}

		for(i=0; i<np+nn; ++i)
			ws[i] /= wsum;

		// grow a tree ...
		grow_rtree(tcodes[ntrees], luts[ntrees], tdepth, tvals, rs, cs, ss, iinds, ws, np+nn, bbox);

		thresholds[ntrees] = -1337.0f;

		++ntrees;

		// update outputs ...
		for(i=0; i<np+nn; ++i)
		{
			float o;

			//
			o = get_tree_output(ntrees-1, rs[i], cs[i], ss[i], iinds[i]);

			//
			os[i] += o;
		}

		// get threshold ...
		threshold = 5.0f;

		do
		{
			//
			threshold -= 0.001f;

			numtps = 0;
			numfps = 0;

			//
			for(i=0; i<np+nn; ++i)
			{
				if( tvals[i]>0 && os[i]>threshold)
					++numtps;
				if(	tvals[i]<0 && os[i]>threshold)
					++numfps;
			}

			//
			tpr = numtps/(float)np;
			fpr = numfps/(float)nn;
		}
		while(tpr<mintpr);

		printf("	** tree %d (%d [s]) ... stage tpr=%f, stage fpr=%f\n", ntrees, (int)(getticks()-t), tpr, fpr);
		fflush(stdout);
	}

	//
	thresholds[ntrees-1] = threshold;

	printf("	** threshold set to %f\n", threshold);

	//
	free(ws);

	//
	return fpr;
}

int find_objects
(
	float rcsq[], int maxndetections,
	int iind,
	float scalefactor, float stridefactor, float minsize, float maxsize
)
{
	float s;
	int ndetections, nrows, ncols;

	//
	nrows = pdims[iind][0];
	ncols = pdims[iind][1];

	ndetections = 0;
	s = minsize;

	while(s<=maxsize)
	{
		float r, c, dr, dc;

		//
		dr = dc = MAX(stridefactor*s, 1.0f);

		//
		for(r=s/2+1; r<=nrows-s/2-1; r+=dr)
			for(c=s/2+1; c<=ncols-s/2-1; c+=dc)
			{
				float q;
				int t;

				t = classify_region(&q, r, c, s, iind);

				if(1==t)
				{
					if(ndetections < maxndetections)
					{
						rcsq[4*ndetections+0] = r;
						rcsq[4*ndetections+1] = c;
						rcsq[4*ndetections+2] = s;
						rcsq[4*ndetections+3] = q;

						//
						++ndetections;
					}
				}
			}

		//
		s = scalefactor*s;
	}

	//
	return ndetections;
}

float get_overlap(float r1, float c1, float s1, float r2, float c2, float s2)
{
	float overr, overc;

	//
	overr = MAX(0, MIN(r1+s1/2, r2+s2/2) - MAX(r1-s1/2, r2-s2/2));
	overc = MAX(0, MIN(c1+s1/2, c2+s2/2) - MAX(c1-s1/2, c2-s2/2));

	//
	return overr*overc/(s1*s1+s2*s2-overr*overc);
}

void search_for_training_data(float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int* np, int* nn, float subsf)
{
	int i, n = 0;

	#define NUMPRNGS 1024
	static int prngsinitialized = 0;
	static uint64_t prngs[NUMPRNGS];

	if(!prngsinitialized)
	{
		// initialize a PRNG for each thread
		for(i=0; i<NUMPRNGS; ++i)
			prngs[i] = 0xFFFF*mwcrand() + 0xFFFF1234FFFF0001LL*mwcrand();

		//
		prngsinitialized = 1;
	}

	*np = 0;
	*nn = 0;

	#pragma omp parallel for
	for(i=0; i<nimages; ++i)
	{
		int thid = omp_get_thread_num();

		#define MAXNDETS 8192
		float dets[4*MAXNDETS];
		int ndets = find_objects(dets, MAXNDETS, i, 1.1f, 0.1f, 24, 1000);

		//printf("%d -> %d %d %d\n", i, ndets, *np, *nn); fflush(stdout);

		int j, k;
		for(j=0; j<ndets; ++j)
		{
			int assigned = 0;

			for(k=contents[i][0]; k<contents[i][1]; ++k)
			{
				float overlap = get_overlap(dets[4*j+0], dets[4*j+1], dets[4*j+2], objects[k][0], objects[k][1], objects[k][2]);

				if(overlap > 0.6f && *np<MAX_N/2 && n<MAX_N)
				#pragma omp critical
				{
					// true positive
					//
					rs[n] = dets[4*j+0];
					cs[n] = dets[4*j+1];
					ss[n] = dets[4*j+2];
					os[n] = dets[4*j+3];

					iinds[n] = i;

					tvals[n] = +1;

					//
					++n;
					++*np;
				}

				if(overlap > 0.4f)
					assigned = 1;
			}

			if(!assigned && (mwcrand_r(&prngs[thid])%1000)/999.0f<subsf && *nn<3*MAX_N/4 && *nn<*np && n<MAX_N)
				#pragma omp critical
				{
					// false positive
					//
					rs[n] = dets[4*j+0];
					cs[n] = dets[4*j+1];
					ss[n] = dets[4*j+2];
					os[n] = dets[4*j+3];

					iinds[n] = i;

					tvals[n] = -1;

					//
					++n;
					++*nn;
				}
		}
	}
}

float sample_training_data(float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int* np, int* nn)
{
	int i, j, n;

	int64_t nw;
	float etpr, efpr;

	int t;

	#define NUMPRNGS 1024
	static int prngsinitialized = 0;
	static uint64_t prngs[NUMPRNGS];

	int stop;

	//
	t = getticks();

	//
	n = 0;

	/*
		object samples
	*/

	if(nobjects > MAX_N)
	{
		printf("* nobjects is too large ... aborting ...\n");
		return -1.0f;
	}

	//
	for(i=0; i<nobjects; ++i)
	{
		//
		int r, c, s, iind;

		iind = objects[i][3];

		r = objects[i][0];
		c = objects[i][1];
		s = objects[i][2];

		//
		if( classify_region(&os[n], r, c, s, iind) == 1 )
		{
			//
			rs[n] = r;
			cs[n] = c;
			ss[n] = s;

			iinds[n] = iind;

			tvals[n] = +1;

			//
			++n;
		}
	}

	*np = n;

	/*
		non-object samples
	*/

	if(!prngsinitialized)
	{
		// initialize a PRNG for each thread
		for(i=0; i<NUMPRNGS; ++i)
			prngs[i] = 0xFFFF*mwcrand() + 0xFFFF1234FFFF0001LL*mwcrand();

		//
		prngsinitialized = 1;
	}

	//
	nw = 0;
	*nn = 0;

	stop = 0;

	#pragma omp parallel
	{
		int thid;

		//
		thid = omp_get_thread_num();

		while(!stop)
		{
			/*
				data mine hard negatives
			*/

			float o;
			int iind, s, r, c;

			// should we sample based on image size?
			iind = mwcrand_r(&prngs[thid])%nimages;

			//if (contents[iind][0] != contents[iind][1])
			//	continue;

			//
			r = mwcrand_r(&prngs[thid])%pdims[iind][0];
			c = mwcrand_r(&prngs[thid])%pdims[iind][1];
			s = objects[mwcrand_r(&prngs[thid])%nobjects][2]; // sample the size of a random object in the pool

			//
			if( classify_region(&o, r, c, s, iind) == 1 )
			{
				// check if the region intersects with a true positive
				// this could probably be done more effciently but we do not expect a large number of objects per image
				int i, ok = 1;
				for(i=contents[iind][0]; i<contents[iind][1]; ++i)
					if(get_overlap(r, c, s, objects[i][0], objects[i][1], objects[i][2]) > 0.5f)
						ok = 0;

				if(ok)
				#pragma omp critical
				{
					//we have a false positive ...
					if(*nn<*np)
					{
						rs[n] = r;
						cs[n] = c;
						ss[n] = s;

						iinds[n] = iind;

						os[n] = o;

						tvals[n] = -1;

						//
						++n;
						++*nn;
					}
					else
						stop = 1;
				}
			}

			if(!stop)
			{
				#pragma omp atomic
				++nw;
			}
		}
	}

	/*
		print the estimated true positive and false positive rates
	*/

	etpr = *np/(float)(nobjects);
	efpr = (float)( *nn/(double)nw );

	printf("* sampling finished ...\n");
	printf("	** elapsed time: %d\n", (int)(getticks()-t));
	printf("	** cascade TPR=%.8f\n", etpr);
	printf("	** cascade FPR=%.8f (%d/%lld)\n", efpr, *nn, (long long int)nw);

	/*
		
	*/

	return efpr;
}

/*
	
*/

static int rs[2*MAX_N];
static int cs[2*MAX_N];
static int ss[2*MAX_N];
static int iinds[2*MAX_N];
static float tvals[2*MAX_N];
static float os[2*MAX_N];

int learn_a_cascade(char* savepath, int8_t bb[], int _tdepth, float stagetpr, float stagefpr, int maxntreesperstage)
{
	int np, nn, i;
	FILE* f = 0;

	//
	bbox[0] = bb[0];
	bbox[1] = bb[1];
	bbox[2] = bb[2];
	bbox[3] = bb[3];

	ntrees = 0;
	tdepth = _tdepth;

	//
	while(1)
	{
		float efpr = sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);

		if(efpr<0.01f)
			break;

		learn_new_stage(stagetpr, stagefpr, maxntreesperstage, tvals, rs, cs, ss, iinds, os, np, nn);

		printf("\n");
	}

	//
	float subsf = 0.01f;

	for(i=0; i<5; ++i)
	{
		//
		printf("* scanning in progress\n");
		search_for_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn, subsf);
		printf("* starting training with np=%d, nn=%d ...\n", np, nn);
		learn_new_stage(stagetpr*stagetpr, stagefpr*stagefpr, maxntreesperstage, tvals, rs, cs, ss, iinds, os, np, nn);
		//
		if(savepath)
		{
			f = fopen(savepath, "wb");
			save_cascade_to_file(f);
			fclose(f);
		}
		// just for estimating FPR for random sampling
		sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
		//
		subsf *= 3;
	}

	//
	printf("* learning process finished\n");
	return 1;
}

/*
	
*/

const char* howto()
{
	return
		"./picolrn <trdata> <cascade-write-path>\n"
	;
}

#ifndef NO_PICOLRN_MAIN
int main(int argc, char* argv[])
{
	FILE* file = 0;

	// initialize the PRNG
	smwcrand(time(0));

	if(argc == 3)
	{
		/*
			* training with default params
			* args: <trdata> <output-cascade>
		*/
		if(!load_training_data(argv[1]))
		{
			printf("* cannot load training data ...\n");
			return 0;
		}

		int8_t bb[] = {-127, +127, -127, +127};
		learn_a_cascade(argv[2], bb, 6, 0.98f, 0.4f, 16);

		file = fopen(argv[2], "wb");
		if(!file || !save_cascade_to_file(file))
		{
			printf("* cannot save result to specified destination\n");
			return 1;
		}
		fclose(file);
	}
	else if(argc == 6)
	{
		/*
			* initializing a new cascade with 0 trees
			* args: <bbox[0]> <bbox[1]> <bbox[2]> <bbox[3]> <tdepth>
		*/
		int tdepth;

		sscanf(argv[1], "%hhd", &bbox[0]);
		sscanf(argv[2], "%hhd", &bbox[1]);
		sscanf(argv[3], "%hhd", &bbox[2]);
		sscanf(argv[4], "%hhd", &bbox[3]);
		sscanf(argv[5], "%d", &tdepth);
		ntrees = 0;

		save_cascade_to_file(stdout);
	}
	else if(argc == 7)
	{
		/*
			* append a new stage to an existing cascade
			* args: <current-cascade> <trdata> <tpr> <fpr> <ntrees> <new-cascade>
		*/
		float stagetpr, stagefpr;
		int maxntreesforstage, np, nn;

		file = fopen(argv[1], "rb");
		if(!file || !load_cascade_from_file(file))
		{
			printf("* cannot load a cascade from '%s'\n", argv[1]);
			return 1;
		}
		fclose(file);

		if(!load_training_data(argv[2]))
		{
			printf("* cannot load the training data from '%s'\n", argv[2]);
			return 1;
		}


		sscanf(argv[3], "%f", &stagetpr);
		sscanf(argv[4], "%f", &stagefpr);
		sscanf(argv[5], "%d", &maxntreesforstage);

		float efpr = sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
		if(efpr < 0.01f)
			search_for_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn, 0.05f);

		printf("* starting training with np=%d, nn=%d ...\n", np, nn);
		learn_new_stage(stagetpr, stagefpr, maxntreesforstage, tvals, rs, cs, ss, iinds, os, np, nn);

		file = fopen(argv[6], "wb");
		if(!file || !save_cascade_to_file(file))
		{
			printf("* cannot save result to specified destination\n");
			return 1;
		}
		fclose(file);
	}
	else
	{
		printf("%s", howto());
	}

	//
	return 0;
}
#endif