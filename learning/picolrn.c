/*
 *	Copyright (c) 2013, Nenad Markus
 *	All rights reserved.
 *
 *	This is an implementation of the algorithm described in the following paper:
 *		N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer,
 *		Object Detection with Pixel Intensity Comparisons Organized in Decision Trees,
 *		http://arxiv.org/abs/1305.4537
 *
 *	Redistribution and use of this program as source code or in binary form, with or without modifications, are permitted provided that the following conditions are met:
 *		1. Redistributions may not be sold, nor may they be used in a commercial product or activity without prior permission from the copyright holder (contact him at nenad.markus@fer.hr).
 *		2. Redistributions may not be used for military purposes.
 *		3. Any published work which utilizes this program shall include the reference to the paper available at http://arxiv.org/abs/1305.4537
 *		4. Redistributions must retain the above copyright notice and the reference to the algorithm on which the implementation is based on, this list of conditions and the following disclaimer.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>

// hyperparameters
#define NRANDS 1024

/*
	auxiliary stuff
*/

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define SQR(x) ((x)*(x))

/*
	- loads an 8-bit grey image saved in the <RID> file format
	- <RID> file contents:
		- a 32-bit signed integer w (image width)
		- a 32-bit signed integer h (image height)
		- an array of w*h unsigned bytes representing pixel intensities
*/

int loadrid(uint8_t* pixels[], int* nrows, int* ncols, const char* path)
{
	FILE* file;
	int w, h;
	
	// open file
	file = fopen(path, "rb");

	if(!file)
	{
		return 0;
	}

	// read width
	fread(&w, sizeof(int), 1, file);
	// read height
	fread(&h, sizeof(int), 1, file);

	// allocate image memory
	*nrows = h;
	*ncols = w;

	*pixels = (uint8_t*)malloc(w*h*sizeof(uint8_t));

	if(!*pixels)
	{
		fclose(file);

		return 0;
	}

	// read image data
	fread(*pixels, sizeof(uint8_t), w*h, file);

	// clean up
	fclose(file);

	// we're done
	return 1;
}

/*
	portable time function
*/

#ifdef __GNUC__
#include <time.h>
float getticks()
{
	struct timespec ts;

	if(clock_gettime(CLOCK_MONOTONIC, &ts) < 0)
	{
		printf("clock_gettime error\n");

		return -1.0f;
	}

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
	regression trees
*/

int bintest(int tcode, int r, int c, int sr, int sc, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	//
	int r1, c1, r2, c2;
	int8_t* p = (int8_t*)&tcode;

	//
	r1 = (256*r + p[0]*sr)/256;
	c1 = (256*c + p[1]*sc)/256;

	r2 = (256*r + p[2]*sr)/256;
	c2 = (256*c + p[3]*sc)/256;

	//
	r1 = MIN(MAX(0, r1), nrows-1);
	c1 = MIN(MAX(0, c1), ncols-1);

	r2 = MIN(MAX(0, r2), nrows-1);
	c2 = MIN(MAX(0, c2), ncols-1);

	//
	return pixels[r1*ldim+c1]<=pixels[r2*ldim+c2];
}

float get_split_error(int tcode, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int i, j;

	double wsum, wsum0, wsum1;
	double wtvalsum0, wtvalsumsqr0, wtvalsum1, wtvalsumsqr1;

	double wmse0, wmse1;

	//
	wsum = wsum0 = wsum1 = wtvalsum0 = wtvalsum1 = wtvalsumsqr0 = wtvalsumsqr1 = 0.0;

	for(i=0; i<indsnum; ++i)
	{
		if( bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
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

int split_training_data(int tcode, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int stop;
	int i, j;

	int n0;

	//
	stop = 0;

	i = 0;
	j = indsnum - 1;

	while(!stop)
	{
		//
		while( !bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
		{
			if( i==j )
				break;
			else
				++i;
		}

		while( bintest(tcode, rs[inds[j]], cs[inds[j]], srs[inds[j]], scs[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]) )
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

	for(i=0; i<indsnum; ++i)
		if( !bintest(tcode, rs[inds[i]], cs[inds[i]], srs[inds[i]], scs[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
			++n0;

	//
	return n0;
}

int grow_subtree(int32_t tcodes[], float lut[], int nodeidx, int d, int maxd, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int inds[], int indsnum)
{
	int i, nrands;

	int32_t tcodelist[2048];
	float spliterrors[2048], bestspliterror;

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

		for(i=0; i<indsnum; ++i)
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
	else if(indsnum <= 1)
	{
		//
		tcodes[nodeidx] = 0;

		//
		grow_subtree(tcodes, lut, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);
		grow_subtree(tcodes, lut, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);

		return 1;
	}

	// generate binary test codes
	nrands = NRANDS;

	for(i=0; i<nrands; ++i)
		tcodelist[i] = mwcrand();

	//
	#pragma omp parallel for
	for(i=0; i<nrands; ++i)
		spliterrors[i] = get_split_error(tcodelist[i], tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);

	//
	bestspliterror = spliterrors[0];
	tcodes[nodeidx] = tcodelist[0];

	for(i=1; i<nrands; ++i)
		if(bestspliterror > spliterrors[i])
		{
			bestspliterror = spliterrors[i];
			tcodes[nodeidx] = tcodelist[i];
		}

	//
	n0 = split_training_data(tcodes[nodeidx], tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, indsnum);

	//
	grow_subtree(tcodes, lut, 2*nodeidx+1, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, &inds[0], n0);
	grow_subtree(tcodes, lut, 2*nodeidx+2, d+1, maxd, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, &inds[n0], indsnum-n0);

	//
	return 1;
}

int grow_rtree(int32_t tcodes[], float lut[], int d, float tvals[], int rs[], int cs[], int srs[], int scs[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], double ws[], int n)
{
	int i;
	int* inds;

	//
	inds = (int*)malloc(n*sizeof(int));

	for(i=0; i<n; ++i)
		inds[i] = i;

	//
	if(!grow_subtree(tcodes, lut, 0, 0, d, tvals, rs, cs, srs, scs, pixelss, nrowss, ncolss, ldims, ws, inds, n))
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
	object samples
*/

#define MAXNUMOS 1000000

static int numos;

static int ors[MAXNUMOS];
static int ocs[MAXNUMOS];
static int oss[MAXNUMOS];

static uint8_t* opixelss[MAXNUMOS];
static int onrowss[MAXNUMOS];
static int oncolss[MAXNUMOS];

int load_object_samples(const char* folder)
{
	char buffer[1024];
	FILE* list;

	//
	printf("Loading object samples from '%s'\n", folder);

	//
	sprintf(buffer, "%s/%s", folder, "list.txt");

	list = fopen(buffer, "r");

	if(!list)
		return 0;

	//
	numos = 0;

	while(fscanf(list, "%s", buffer) == 1) // read image name
	{
		char fullpath[1024];

		int nrows, ncols;
		uint8_t* opixels;
		int i, n;

		//
		if(numos >= MAXNUMOS)
		{
			printf("maximum allowed number of object samples exceeded: terminating ...\n");

			return 0;
		}

		// load image
		sprintf(fullpath, "%s/%s", folder, buffer);

		if(!loadrid(&opixels, &nrows, &ncols, fullpath))
			return 0;

		// number of samples associated with this image
		if(fscanf(list, "%d", &n) != 1)
			return 0;

		// get samples
		for(i=0; i<n; ++i)
		{
			int r, c, s;

			//
			if(fscanf(list, "%d %d %d", &r, &c, &s) != 3)
				return 0;

			//
			ors[numos] = r;
			ocs[numos] = c;
			oss[numos] = s;

			opixelss[numos] = opixels;
			onrowss[numos] = nrows;
			oncolss[numos] = ncols;

			//
			++numos;
		}
	}

	fclose(list);

	return 1;
}

/*
	background images
*/

#define MAXNUMBS 100000

static int numbs = 0;

static uint8_t* bpixelss[MAXNUMBS];
int bnrowss[MAXNUMBS];
int bncolss[MAXNUMBS];

int load_background_images(char* folder)
{
	FILE* list;
	char path[1024], name[1024];

	//
	printf("Loading background images from '%s'\n", folder);

	//
	sprintf(path, "%s/list.txt", folder);

	list = fopen(path, "r");

	if(!list)
		return 0;

	//
	while(fscanf(list, "%s", name)==1 && numbs<MAXNUMBS)
	{
		//
		sprintf(path, "%s/%s", folder, name);

		//
		if(loadrid(&bpixelss[numbs], &bnrowss[numbs], &bncolss[numbs], path))
			++numbs;
	}

	//
	fclose(list);

	//
	return 1;
}

/*
	
*/

float tsr, tsc;
int tdepth;
int ntrees=0;

int32_t tcodes[1024][1024];
float luts[1024][1024];

float thresholds[1024];

/*
	
*/

int load_from_file(const char* path)
{
	int i;
	FILE* file;

	//
	file = fopen(path, "rb");

	if(!file)
		return 0;

	//
	fread(&tsr, sizeof(float), 1, file);
	fread(&tsc, sizeof(float), 1, file);

	fread(&tdepth, sizeof(int), 1, file);

	fread(&ntrees, sizeof(int), 1, file);

	//
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

int save_to_file(const char* path)
{
	int i;
	FILE* file;

	//
	file = fopen(path, "wb");

	if(!file)
		return 0;

	//
	fwrite(&tsr, sizeof(float), 1, file);
	fwrite(&tsc, sizeof(float), 1, file);

	fwrite(&tdepth, sizeof(int), 1, file);

	fwrite(&ntrees, sizeof(int), 1, file);

	//
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

float get_tree_output(int i, int r, int c, int sr, int sc, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int idx, j;

	//
	idx = 1;

	for(j=0; j<tdepth; ++j)
		idx = 2*idx + bintest(tcodes[i][idx-1], r, c, sr, sc, pixels, nrows, ncols, ldim);

	//
	return luts[i][idx - (1<<tdepth)];
}

int classify_region(float* o, int r, int c, int s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int i, j, sr, sc;

	//
	if(!ntrees)
		return 1;

	//
	sr = (int)(tsr*s);
	sc = (int)(tsc*s);

	*o = 0.0f;

	//
	for(i=0; i<ntrees; ++i)
	{
		//
		*o += get_tree_output(i, r, c, sr, sc, pixels, nrows, ncols, ldim);

		//
		if(*o <= thresholds[i])
			return -1;

		//
		++i;
	}

	//
	return 1;
}

int learn_new_stage(float mintpr, float maxfpr, int maxntrees, int classs[], int rs[], int cs[], int ss[], uint8_t* ppixels[], int nrowss[], int ncolss[], float os[], int np, int nn)
{
	int i;

	float* tvals;

	int* srs;
	int* scs;

	double* ws;
	double wsum;

	float threshold, tpr, fpr;

	//
	tvals = (float*)malloc((np+nn)*sizeof(float));

	for(i=0; i<np+nn; ++i)
		if(classs[i])
			tvals[i] = +1.0f; // object
		else
			tvals[i] = -1.0f; // non-object

	//
	srs = (int*)malloc((np+nn)*sizeof(int));
	scs = (int*)malloc((np+nn)*sizeof(int));

	for(i=0; i<np+nn; ++i)
	{
		srs[i] = (int)( tsr*ss[i] );
		scs[i] = (int)( tsc*ss[i] );
	}

	//
	ws = (double*)malloc((np+nn)*sizeof(double));

	//
	maxntrees = ntrees + maxntrees;
	fpr = 1.0f;

	while(ntrees<maxntrees && fpr>maxfpr)
	{
		float t;
		int numtps, numfps;

		// compute weights ...
		wsum = 0.0;

		for(i=0; i<np+nn; ++i)
		{
			if(classs[i])
				ws[i] = exp(-1.0*os[i])/np;
			else
				ws[i] = exp(+1.0*os[i])/nn;

			wsum += ws[i];
		}

		for(i=0; i<np+nn; ++i)
			ws[i] /= wsum;

		// grow a tree ...
		t = getticks();

		grow_rtree(tcodes[ntrees], luts[ntrees], tdepth, tvals, rs, cs, srs, scs, ppixels, nrowss, ncolss, ncolss, ws, np+nn);

		printf("\r");
		printf("	tree %d (%f [s]) ...", ntrees+1, getticks()-t);

		thresholds[ntrees] = -1337.0f;

		++ntrees;

		// update outputs ...
		for(i=0; i<np+nn; ++i)
		{
			float o;

			//
			o = get_tree_output(ntrees-1, rs[i], cs[i], srs[i], scs[i], ppixels[i], nrowss[i], ncolss[i], ncolss[i]);

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
				if( classs[i] && os[i]>threshold)
					++numtps;
				if(!classs[i] && os[i]>threshold)
					++numfps;
			}

			//
			tpr = numtps/(float)np;
			fpr = numfps/(float)nn;
		}
		while(tpr<mintpr);

		printf(" tpr=%f, fpr=%f\t", tpr, fpr);
		fflush(stdout);
	}

	//
	thresholds[ntrees-1] = threshold;

	printf("\n");
	printf("	threshold set to %f\n", threshold);

	//
	free(tvals);

	free(srs);
	free(scs);

	free(ws);

	//
	return 1;
}

float sample_training_data(int classs[], int rs[], int cs[], int ss[], uint8_t* ppixels[], int nrowss[], int ncolss[], float os[], int maxn, int* np, int* nn)
{
	int i, n;

	int64_t nw;
	float etpr, efpr;

	int t;

	#define NUMPRNGS 1024
	static int prngsinitialized = 0;
	static int64_t prngs[NUMPRNGS];

	int stop;

	//
	printf("- sampling training data (randomized) ...\n");

	t = getticks();

	//
	n = 0;

	/*
		object samples
	*/

	for(i=0; i<numos; ++i)
		if( classify_region(&os[n], ors[i], ocs[i], oss[i], opixelss[i], onrowss[i], oncolss[i], oncolss[i]) == 1 )
		{
			//
			rs[n] = ors[i];
			cs[n] = ocs[i];
			ss[n] = oss[i];

			ppixels[n] = opixelss[i];
			nrowss[n] = onrowss[i];
			ncolss[n] = oncolss[i];

			classs[n] = +1;

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

	#pragma omp parallel
	{
		int thid;

		//
		thid = omp_get_thread_num();

		while(!stop)
		{
			float o;
			int idx, s, r, c, nrows, ncols;
			uint8_t* pixels;

			//
			idx = mwcrand_r(&prngs[thid])%numbs;

			//
			pixels = bpixelss[idx];
			nrows = bnrowss[idx];
			ncols = bncolss[idx];

			r = mwcrand_r(&prngs[thid])%bnrowss[idx];
			c = mwcrand_r(&prngs[thid])%bncolss[idx];
			s = mwcrand_r(&prngs[thid])%( 2*MIN(MIN(r, nrows-r), MIN(c, ncols-c)) + 1 );

			if(s<24)
				continue;

			//
			if( classify_region(&o, r, c, s, pixels, nrows, ncols, ncols) == 1 )
			{
				//we have a false positive ...
				#pragma omp critical
				{
					if(n<maxn)
					{
						rs[n] = r;
						cs[n] = c;
						ss[n] = s;

						ppixels[n] = pixels;
						nrowss[n] = nrows;
						ncolss[n] = ncols;

						os[n] = o;

						classs[n] = 0;

						//
						++n;
						++*nn;
					}
					else
						stop = 1;
				}
			}

			#pragma omp atomic
			++nw;
		}
	}

	/*
		print estimated true positive and false positive rates
	*/

	etpr = *np/(float)numos;
	efpr = (float)( *nn/(double)nw );

	printf("	tpr (sampling): %.8f\n", etpr);
	printf("	fpr (sampling): %.8f (%d/%lld)\n", efpr, *nn, (long long int)nw);
	printf("	elapsed time: %f  [sec]\n", getticks()-t);

	/*
		
	*/

	return efpr;
}

int append_stages(char* src, char* dst, int maxnumstagestoappend, float targetfpr, float minstagetpr, float maxstagefpr, int maxnumtreesperstage)
{
	#define MAXMAXNUMSAMPLES 2*MAXNUMOS

	static int rs[MAXMAXNUMSAMPLES];
	static int cs[MAXMAXNUMSAMPLES];
	static int ss[MAXMAXNUMSAMPLES];

	static int classs[MAXMAXNUMSAMPLES];

	static uint8_t* ppixels[MAXMAXNUMSAMPLES];
	static int nrowss[MAXMAXNUMSAMPLES];
	static int ncolss[MAXMAXNUMSAMPLES];

	static float os[MAXMAXNUMSAMPLES];

	//
	int i, maxnumsamples, np, nn;

	//
	if(!load_from_file(src))
		return 0;

	//
	maxnumsamples = 2*numos;

	np = nn = 0;

	for(i=0; i<maxnumstagestoappend; ++i)
	{
		float currentfpr;

		printf("--------------------------------------------------------------------\n");
		printf("%d/%d\n", i+1, maxnumstagestoappend);

		/*
			sample training set
		*/

		printf("\n");

		currentfpr = sample_training_data(classs, rs, cs, ss, ppixels, nrowss, ncolss, os, maxnumsamples, &np, &nn);

		if(currentfpr <= targetfpr)
		{
			printf("- target FPR achieved ... terminating the learning process ...\n");

			break;
		}

		/*
			learn decision trees in the current stage
		*/

		printf("\n");
		printf("- learning a new stage ...\n");
		printf("	npositives: %d, nnegatives: %d\n", np, nn);

		learn_new_stage(minstagetpr, maxstagefpr, maxnumtreesperstage, classs, rs, cs, ss, ppixels, nrowss, ncolss, os, np, nn);

		/*
			
		*/

		//
		printf("\n");

		if(save_to_file(dst))
			printf("- saving partial results to '%s' ...\n", dst);
		else
			printf("- saving results to '%s' has failed ...\n", dst);

		printf("\n");
	}

	printf("--------------------------------------------------------------------\n");

	//
	return 1;
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
	float t;

	char* src;
	char* dst;
	char* objspath;
	char* nonobjimgspath;

	int maxnstages;
	float targetfpr;
	float minstagetpr, maxstagefpr;
	int tdepths, maxnumtreesperstage;

	//
	if(argc == 5)
	{
		sscanf(argv[1], "%f", &tsr);
		sscanf(argv[2], "%f", &tsc);

		sscanf(argv[3], "%d", &tdepth);

		//
		ntrees = 0;

		//
		if(!save_to_file(argv[4]))
			return 0;

		//
		printf("INITIALIZING: (%f, %f, %d)\n", tsr, tsc, tdepth);

		//
		return 0;
	}
	else if(argc == 10)
	{
		src = argv[1];

		objspath = argv[2];
		nonobjimgspath = argv[3];

		sscanf(argv[4], "%d", &maxnstages);
		sscanf(argv[5], "%f", &targetfpr);
		sscanf(argv[6], "%f", &minstagetpr);
		sscanf(argv[7], "%f", &maxstagefpr);
		sscanf(argv[8], "%d", &maxnumtreesperstage);

		dst = argv[9];
	}
	else
	{
		printf("%s", howto());
		return 0;
	}

	// initialize PRNG
	smwcrand(time(0));

	//
	t = getticks();
	if(!load_object_samples(objspath))
	{
		printf("cannot load object samples ... exiting ...\n");
		return 1;
	}
	printf("%d object samples loaded in %f [s]\n", numos, getticks()-t);

	//
	t = getticks();
	if(!load_background_images(nonobjimgspath))
	{
		printf("cannot load background images ... exiting ...\n");
		return 1;
	}
	printf("%d background images loaded in %f [s]\n", numbs, getticks()-t);

	//
	t = getticks();
	printf("LEARNING ...\n");
	append_stages(src, dst, maxnstages, targetfpr, minstagetpr, maxstagefpr, maxnumtreesperstage);
	printf("FINISHED ...\n");

	printf("elapsed time: %f [sec]\n", getticks()-t);

	printf("\n");

	//
	return 0;
}
