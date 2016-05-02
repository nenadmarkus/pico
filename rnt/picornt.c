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

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))

#include <stdint.h>

/*
	
*/

int run_cascade(void* cascade, float* o, int r, int c, int s, void* vppixels, int nrows, int ncols, int ldim)
{
	//
	int i, j, idx;

	uint8_t* pixels;

	int tdepth, ntrees, offset;

	int8_t* ptree;
	int8_t* tcodes;
	float* lut;
	float thr;

	//
	pixels = (uint8_t*)vppixels;

	//
	tdepth = ((int*)cascade)[2];
	ntrees = ((int*)cascade)[3];

	//
	r = r*256;
	c = c*256;

	if( (r+128*s)/256>=nrows || (r-128*s)/256<0 || (c+128*s)/256>=ncols || (c-128*s)/256<0 )
		return -1;

	//
	offset = ((1<<tdepth)-1)*sizeof(int32_t) + (1<<tdepth)*sizeof(float) + 1*sizeof(float);
	ptree = (int8_t*)cascade + 2*sizeof(float) + 2*sizeof(int);

	*o = 0.0f;

	for(i=0; i<ntrees; ++i)
	{
		//
		tcodes = ptree - 4;
		lut = (float*)(ptree + ((1<<tdepth)-1)*sizeof(int32_t));
		thr = *(float*)(ptree + ((1<<tdepth)-1)*sizeof(int32_t) + (1<<tdepth)*sizeof(float));

		//
		idx = 1;

		for(j=0; j<tdepth; ++j)
			idx = 2*idx + (pixels[(r+tcodes[4*idx+0]*s)/256*ldim+(c+tcodes[4*idx+1]*s)/256]<=pixels[(r+tcodes[4*idx+2]*s)/256*ldim+(c+tcodes[4*idx+3]*s)/256]);

		*o = *o + lut[idx-(1<<tdepth)];

		//
		if(*o<=thr)
			return -1;
		else
			ptree = ptree + offset;
	}

	//
	*o = *o - thr;

	return +1;
}

int run_rotated_cascade(void* cascade, float* o, int r, int c, int s, float a, void* vppixels, int nrows, int ncols, int ldim)
{
	//
	int i, j, idx;

	uint8_t* pixels;

	int tdepth, ntrees, offset;

	int8_t* ptree;
	int8_t* tcodes;
	float* lut;
	float thr;

	static int qcostable[32+1] = {256, 251, 236, 212, 181, 142, 97, 49, 0, -49, -97, -142, -181, -212, -236, -251, -256, -251, -236, -212, -181, -142, -97, -49, 0, 49, 97, 142, 181, 212, 236, 251, 256};
	static int qsintable[32+1] = {0, 49, 97, 142, 181, 212, 236, 251, 256, 251, 236, 212, 181, 142, 97, 49, 0, -49, -97, -142, -181, -212, -236, -251, -256, -251, -236, -212, -181, -142, -97, -49, 0};

	//
	pixels = (uint8_t*)vppixels;

	//
	tdepth = ((int*)cascade)[2];
	ntrees = ((int*)cascade)[3];

	//
	r = r*65536;
	c = c*65536;

	if( (r+46341*s)/65536>=nrows || (r-46341*s)/65536<0 || (c+46341*s)/65536>=ncols || (c-46341*s)/65536<0 )
		return -1;

	//
	offset = ((1<<tdepth)-1)*sizeof(int32_t) + (1<<tdepth)*sizeof(float) + 1*sizeof(float);
	ptree = (int8_t*)cascade + 2*sizeof(float) + 2*sizeof(int);

	*o = 0.0f;

	int qsin = s*qsintable[(int)(32*a)]; //s*(int)(256.0f*sinf(2*M_PI*a));
	int qcos = s*qcostable[(int)(32*a)]; //s*(int)(256.0f*cosf(2*M_PI*a));

	for(i=0; i<ntrees; ++i)
	{
		//
		tcodes = ptree - 4;
		lut = (float*)(ptree + ((1<<tdepth)-1)*sizeof(int32_t));
		thr = *(float*)(ptree + ((1<<tdepth)-1)*sizeof(int32_t) + (1<<tdepth)*sizeof(float));

		//
		idx = 1;

		for(j=0; j<tdepth; ++j)
		{
			int r1, c1, r2, c2;

			//
			r1 = (r + qcos*tcodes[4*idx+0] - qsin*tcodes[4*idx+1])/65536;
			c1 = (c + qsin*tcodes[4*idx+0] + qcos*tcodes[4*idx+1])/65536;

			r2 = (r + qcos*tcodes[4*idx+2] - qsin*tcodes[4*idx+3])/65536;
			c2 = (c + qsin*tcodes[4*idx+2] + qcos*tcodes[4*idx+3])/65536;

			//
			idx = 2*idx + (pixels[r1*ldim+c1]<=pixels[r2*ldim+c2]);
		}

		*o = *o + lut[idx-(1<<tdepth)];

		//
		if(*o<=thr)
			return -1;
		else
			ptree = ptree + offset;
	}

	//
	*o = *o - thr;

	return +1;
}

int find_objects
		(
			float rs[], float cs[], float ss[], float qs[], int maxndetections,
			void* cascade, float angle, // * `angle` is a number between 0 and 1 that determines the counterclockwise in-plane rotation of the cascade: 0.0f corresponds to 0 radians and 1.0f corresponds to 2*pi radians
			void* pixels, int nrows, int ncols, int ldim,
			float scalefactor, float stridefactor, float minsize, float maxsize
		)
{
	float s;
	int ndetections;

	//
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

				if(0.0f==angle)
					t = run_cascade(cascade, &q, r, c, s, pixels, nrows, ncols, ldim);
				else
					t = run_rotated_cascade(cascade, &q, r, c, s, angle, pixels, nrows, ncols, ldim);

				if(1==t)
				{
					if(ndetections < maxndetections)
					{
						qs[ndetections] = q;
						rs[ndetections] = r;
						cs[ndetections] = c;
						ss[ndetections] = s;

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

/*
	
*/

float get_overlap(float r1, float c1, float s1, float r2, float c2, float s2)
{
	float overr, overc;

	//
	overr = MAX(0, MIN(r1+s1/2, r2+s2/2) - MAX(r1-s1/2, r2-s2/2));
	overc = MAX(0, MIN(c1+s1/2, c2+s2/2) - MAX(c1-s1/2, c2-s2/2));

	//
	return overr*overc/(s1*s1+s2*s2-overr*overc);
}

void ccdfs(int a[], int i, float rs[], float cs[], float ss[], int n)
{
	int j;

	//
	for(j=0; j<n; ++j)
		if(a[j]==0 && get_overlap(rs[i], cs[i], ss[i], rs[j], cs[j], ss[j])>0.3f)
		{
			//
			a[j] = a[i];

			//
			ccdfs(a, j, rs, cs, ss, n);
		}
}

int find_connected_components(int a[], float rs[], float cs[], float ss[], int n)
{
	int i, ncc, cc;

	//
	if(!n)
		return 0;

	//
	for(i=0; i<n; ++i)
		a[i] = 0;

	//
	ncc = 0;
	cc = 1;

	for(i=0; i<n; ++i)
		if(a[i] == 0)
		{
			//
			a[i] = cc;

			//
			ccdfs(a, i, rs, cs, ss, n);

			//
			++ncc;
			++cc;
		}

	//
	return ncc;
}

int cluster_detections(float rs[], float cs[], float ss[], float qs[], int n)
{
	int idx, ncc, cc;
	int a[4096];

	//
	ncc = find_connected_components(a, rs, cs, ss, n);

	if(!ncc)
		return 0;

	//
	idx = 0;

	for(cc=1; cc<=ncc; ++cc)
	{
		int i, k;

		float sumqs=0.0f, sumrs=0.0f, sumcs=0.0f, sumss=0.0f;

		//
		k = 0;

		for(i=0; i<n; ++i)
			if(a[i] == cc)
			{
				sumqs += qs[i];
				sumrs += rs[i];
				sumcs += cs[i];
				sumss += ss[i];

				++k;
			}

		//
		qs[idx] = sumqs; // accumulated confidence measure

		//
		rs[idx] = sumrs/k;
		cs[idx] = sumcs/k;
		ss[idx] = sumss/k;

		//
		++idx;
	}

	//
	return idx;
}