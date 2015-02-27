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
#include <stdint.h>
#include <math.h>

#define MAX(a, b) ((a)>(b)?(a):(b))
#define ABS(x) ((x)>0?(x):(-(x)))

/*
	
*/

float tsr, tsc;
int tdepth;
int ntrees=0;

int32_t tcodes[4096][1024];
float luts[4096][1024];

float thresholds[4096];

/*
	
*/

int load_cascade(const char* path)
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

int save_cascade(const char* path)
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

void print_c_code(const char* name, float rotation)
{
	int i, j, maxr, maxc;

	int qsin, qcos, q;

	static int16_t rtcodes[4096][1024][4];

	// generate rotated binary tests
	q = (1<<16);

	qsin = (int)( q*sin(rotation) );
	qcos = (int)( q*cos(rotation) );

	maxr = 0;
	maxc = 0;

	for(i=0; i<ntrees; ++i)
		for(j=0; j<(1<<tdepth)-1; ++j)
		{
			int8_t* p = (int8_t*)&tcodes[i][j];

			//
			rtcodes[i][j][0] = (p[0]*qcos - p[1]*qsin)/q;
			rtcodes[i][j][1] = (p[0]*qsin + p[1]*qcos)/q;

			rtcodes[i][j][2] = (p[2]*qcos - p[3]*qsin)/q;
			rtcodes[i][j][3] = (p[2]*qsin + p[3]*qcos)/q;

			//
			maxr = MAX(maxr, MAX(ABS(rtcodes[i][j][0]), ABS(rtcodes[i][j][2])));
			maxc = MAX(maxc, MAX(ABS(rtcodes[i][j][1]), ABS(rtcodes[i][j][3])));
		}

	//
	printf("int %s(float* o, int r, int c, int s, void* vppixels, int nrows, int ncols, int ldim)\n", name);
	printf("{\n");

	//
	printf("	int i, idx, sr, sc;\n");
	printf("	uint8_t* pixels;\n");

	//
	printf("\n");
	printf("	static int16_t tcodes[%d][%d][4] =\n", ntrees, 1<<tdepth);
	printf("	{\n");
	for(i=0; i<ntrees; ++i)
	{
		printf("		{{0, 0, 0, 0}");
		for(j=0; j<(1<<tdepth)-1; ++j)
			printf(", {%d, %d, %d, %d}", rtcodes[i][j][0], rtcodes[i][j][1], rtcodes[i][j][2], rtcodes[i][j][3]);
		printf("},\n");
	}
	printf("	};\n");

	//
	printf("\n");
	printf("	static float lut[%d][%d] =\n", ntrees, 1<<tdepth);
	printf("	{\n");
	for(i=0; i<ntrees; ++i)
	{
		printf("		{");
		for(j=0; j<(1<<tdepth)-1; ++j)
			printf("%ff, ", luts[i][j]);
		printf("%ff},\n", luts[i][(1<<tdepth)-1]);
	}
	printf("	};\n");

	//
	printf("\n");
	printf("	static float thresholds[%d] =\n", ntrees);
	printf("	{\n\t\t");
	for(i=0; i<ntrees-1; ++i)
		printf("%ff, ", thresholds[i]);
	printf("%ff\n", thresholds[ntrees-1]);
	printf("	};\n");

	//
	printf("\n");
	printf("	sr = (int)(%ff*s);\n", tsr);
	printf("	sc = (int)(%ff*s);\n", tsc);

	// generate the code that checks image boundaries
	printf("\n");
	printf("	if( (256*r+%d*sr)/256>=nrows || (256*r-%d*sr)/256<0 || (256*c+%d*sc)/256>=ncols || (256*c-%d*sc)/256<0 )\n", maxr, maxr, maxc, maxc);
	printf("		return -1;\n");

	printf("\n");
	printf("	pixels = (uint8_t*)vppixels;\n");

	//
	printf("\n");
	printf("	*o = 0.0f;\n\n");
	///printf("	pixels = &pixels[r*ldim+c];\n");
	printf("	for(i=0; i<%d; ++i)\n", ntrees);
	printf("	{\n");
	printf("		idx = 1;\n");
	for(i=0; i<tdepth; ++i)
	{
		printf("		idx = 2*idx + (pixels[(256*r+tcodes[i][idx][0]*sr)/256*ldim + (256*c+tcodes[i][idx][1]*sc)/256]<=pixels[(256*r+tcodes[i][idx][2]*sr)/256*ldim + (256*c+tcodes[i][idx][3]*sc)/256]);\n");
		///printf("		idx = 2*idx + (pixels[tcodes[i][idx][0]*sr/256*ldim + tcodes[i][idx][1]*sc/256]<=pixels[tcodes[i][idx][2]*sr/256*ldim + tcodes[i][idx][3]*sc/256]);\n");
	}
	printf("\n		*o = *o + lut[i][idx-%d];\n\n", 1<<tdepth);
	printf("		if(*o<=thresholds[i])\n\t\t\treturn -1;\n");
	printf("	}\n");

	printf("\n	*o = *o - thresholds[%d];\n", ntrees-1);
	printf("\n");
	printf("	return +1;\n");

	//
	printf("}\n");
}

/*
	
*/

int main(int argc, char* argv[])
{
	//
	if(argc == 3)
	{
		load_cascade(argv[1]);
		print_c_code(argv[2], 0.0f);
	}
	else if(argc == 4)
	{
		float rotation;

		load_cascade(argv[1]);
		sscanf(argv[2], "%f", &rotation);
		print_c_code(argv[3], rotation);
	}
	else
	{
		printf("* specify arguments: <cascade> <in-plane rotation> <detection function name>\n");
		return 0;
	}

	//
	return 0;
}
