/*
 *  This code is released under the MIT License.
 *  Copyright (c) 2013 Nenad Markus
 */

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>

#include <omp.h>

// hyperparameters
#define NRANDS 1024

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

#define MAX_N 2000000

int N = 0;
uint8_t* ppixels[MAX_N];
int pdims[MAX_N][2]; // (nrows, ncols)

int nbackground = 0;
int background[MAX_N]; // i

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
	N = 0;

	nbackground = 0;
	nobjects = 0;

	while( load_image(&ppixels[N], &pdims[N][0], &pdims[N][1], file) )
	{
		int i, n;

		//
		if(fread(&n, sizeof(int), 1, file) != 1)
			return 1;

		if(!n)
		{
			background[nbackground] = N;
			++nbackground;
		}
		else
		{
			for(i=0; i<n; ++i)
			{
				fread(&objects[nobjects][0], sizeof(int), 1, file); // r
				fread(&objects[nobjects][1], sizeof(int), 1, file); // c
				fread(&objects[nobjects][2], sizeof(int), 1, file); // s

				objects[nobjects][3] = N; // i

				//
				++nobjects;
			}
		}

		//
		++N;
	}

	//
	return 1;
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
	int i, j;

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

int load_cascade_from_file(const char* path)
{
	int i;
	FILE* file;

	//
	file = fopen(path, "rb");

	if(!file)
		return 0;

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
	fclose(file);

	//
	return 1;
}

int save_cascade_to_file(const char* path)
{
	int i;
	FILE* file;

	//
	file = fopen(path, "wb");

	if(!file)
		return 0;

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
	fclose(file);

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
	int i, sr, sc;

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

int learn_new_stage(float mintpr, float maxfpr, int maxntrees, float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int np, int nn)
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
			threshold -= 0.005f;

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
	return 1;
}

float sample_training_data(float tvals[], int rs[], int cs[], int ss[], int iinds[], float os[], int* np, int* nn)
{
	int i, n;

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

	for(i=0; i<nobjects; ++i)
		if( classify_region(&os[n], objects[i][0], objects[i][1], objects[i][2], objects[i][3]) == 1 )
		{
			//
			rs[n] = objects[i][0];
			cs[n] = objects[i][1];
			ss[n] = objects[i][2];


			iinds[n] = objects[i][3];

			tvals[n] = +1;

			//
			++n;
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

	if(nbackground)
	{
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
				int iind, s, r, c, nrows, ncols;
				uint8_t* pixels;

				//
				iind = background[ mwcrand_r(&prngs[thid])%nbackground ];

				//
				r = mwcrand_r(&prngs[thid])%pdims[iind][0];
				c = mwcrand_r(&prngs[thid])%pdims[iind][1];
				s = objects[mwcrand_r(&prngs[thid])%nobjects][2]; // sample the size of a random object in the pool

				//
				if( classify_region(&o, r, c, s, iind) == 1 )
				{
					//we have a false positive ...
					#pragma omp critical
					{
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
	}
	else
		nw = 1;

	/*
		print the estimated true positive and false positive rates
	*/

	etpr = *np/(float)nobjects;
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

int learn_with_default_parameters(char* trdata, char* dst)
{
	int i, np, nn;
	float fpr;

	//
	if(!load_training_data(trdata))
	{
		printf("* cannot load training data ...\n");
		return 0;
	}

	//
	bbox[0] = -127;
	bbox[1] = +127;
	bbox[2] = -127;
	bbox[3] = +127;

	tdepth = 5;

	if(!save_cascade_to_file(dst))
			return 0;

	//
	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
	learn_new_stage(0.9800f, 0.5f, 4, tvals, rs, cs, ss, iinds, os, np, nn);
	save_cascade_to_file(dst);

	printf("\n");

	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
	learn_new_stage(0.9850f, 0.5f, 8, tvals, rs, cs, ss, iinds, os, np, nn);
	save_cascade_to_file(dst);

	printf("\n");

	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
	learn_new_stage(0.9900f, 0.5f, 16, tvals, rs, cs, ss, iinds, os, np, nn);
	save_cascade_to_file(dst);

	printf("\n");

	sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
	learn_new_stage(0.9950f, 0.5f, 32, tvals, rs, cs, ss, iinds, os, np, nn);
	save_cascade_to_file(dst);

	printf("\n");

	//
	while(sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn) > 1e-6f)
	{
		learn_new_stage(0.9975f, 0.5f, 64, tvals, rs, cs, ss, iinds, os, np, nn);
		save_cascade_to_file(dst);

		printf("\n");
	}

	//
	printf("* target FPR achieved ... terminating the learning process ...\n");
}

/*
	
*/

const char* howto()
{
	return
		"TODO\n"
	;
}

int main(int argc, char* argv[])
{
	// initialize the PRNG
	smwcrand(time(0));

	//
	if(argc == 3)
	{
		learn_with_default_parameters(argv[1], argv[2]);
	}
	else if(argc == 7)
	{
		int dummy;

		//
		sscanf(argv[1], "%d", &dummy); bbox[0] = dummy;
		sscanf(argv[2], "%d", &dummy); bbox[1] = dummy;
		sscanf(argv[3], "%d", &dummy); bbox[2] = dummy;
		sscanf(argv[4], "%d", &dummy); bbox[3] = dummy;

		//
		sscanf(argv[5], "%d", &tdepth);

		//
		ntrees = 0;

		//
		if(!save_cascade_to_file(argv[6]))
			return 0;

		//
		printf("* initializing:\n");
		printf("	** bbox = (%d, %d, %d, %d)\n", bbox[0], bbox[1], bbox[2], bbox[3]);
		printf("	** tdepth = %d\n", tdepth);

		//
		return 0;
	}
	else if(argc == 7)
	{
		float tpr, fpr;
		int ntrees, np, nn;

		//
		if(!load_cascade_from_file(argv[1]))
		{
			printf("* cannot load a cascade from '%s'\n", argv[1]);
			return 1;
		}

		if(!load_training_data(argv[2]))
		{
			printf("* cannot load the training data from '%s'\n", argv[2]);
			return 1;
		}

		//
		sscanf(argv[3], "%f", &tpr);
		sscanf(argv[4], "%f", &fpr);
		sscanf(argv[5], "%d", &ntrees);

		//
		sample_training_data(tvals, rs, cs, ss, iinds, os, &np, &nn);
		learn_new_stage(tpr, fpr, ntrees, tvals, rs, cs, ss, iinds, os, np, nn);

		//
		if(!save_cascade_to_file(argv[6]))
			return 1;
	}
	else
	{
		printf("%s", howto());
		return 0;
	}

	//
	return 0;
}
