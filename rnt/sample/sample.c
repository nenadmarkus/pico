/*
 *	Copyright (c) 2013, Nenad Markus
 *	All rights reserved.
 *
 *	This is an implementation of the algorithm described in the following paper:
 *		N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer,
 *		A method for object detection based on pixel intensity comparisons,
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

#include <cv.h>
#include <highgui.h>

//
#include "../picornt.h"

/*
	a portable time function
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
	
*/

void* cascade = 0;

int minsize;
int maxsize;

float scalefactor;
float stridefactor;

float qthreshold;

int usepyr;
int noclustering;
int verbose;

void process_image(IplImage* frame, int draw)
{
	int i, j;
	float t;

	uint8_t* pixels;
	int nrows, ncols, ldim;

	#define MAXNDETECTIONS 2048
	int ndetections;
	float qs[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];

	static IplImage* gray = 0;
	static IplImage* pyr[5] = {0, 0, 0, 0, 0};

	/*
		...
	*/

	//
	if(!pyr[0])
	{
		//
		gray = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, 1);

		//
		pyr[0] = gray;
		pyr[1] = cvCreateImage(cvSize(frame->width/2, frame->height/2), frame->depth, 1);
		pyr[2] = cvCreateImage(cvSize(frame->width/4, frame->height/4), frame->depth, 1);
		pyr[3] = cvCreateImage(cvSize(frame->width/8, frame->height/8), frame->depth, 1);
		pyr[4] = cvCreateImage(cvSize(frame->width/16, frame->height/16), frame->depth, 1);
	}

	// get grayscale image
	if(frame->nChannels == 3)
		cvCvtColor(frame, gray, CV_RGB2GRAY);
	else
		cvCopy(frame, gray, 0);

	// perform detection with the pico library
	t = getticks();

	if(usepyr)
	{
		int nd;

		//
		pyr[0] = gray;

		pixels = (uint8_t*)pyr[0]->imageData;
		nrows = pyr[0]->height;
		ncols = pyr[0]->width;
		ldim = pyr[0]->widthStep;

		ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, cascade, 0.0f, pixels, nrows, ncols, ldim, scalefactor, stridefactor, MAX(16, minsize), MIN(128, maxsize));

		for(i=1; i<5; ++i)
		{
			cvResize(pyr[i-1], pyr[i], CV_INTER_LINEAR);

			pixels = (uint8_t*)pyr[i]->imageData;
			nrows = pyr[i]->height;
			ncols = pyr[i]->width;
			ldim = pyr[i]->widthStep;

			nd = find_objects(&rs[ndetections], &cs[ndetections], &ss[ndetections], &qs[ndetections], MAXNDETECTIONS-ndetections, cascade, 0.0f, pixels, nrows, ncols, ldim, scalefactor, stridefactor, MAX(64, minsize>>i), MIN(128, maxsize>>i));

			for(j=ndetections; j<ndetections+nd; ++j)
			{
				rs[j] = (1<<i)*rs[j];
				cs[j] = (1<<i)*cs[j];
				ss[j] = (1<<i)*ss[j];
			}

			ndetections = ndetections + nd;
		}
	}
	else
	{
		//
		pixels = (uint8_t*)gray->imageData;
		nrows = gray->height;
		ncols = gray->width;
		ldim = gray->widthStep;

		//
		ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, cascade, 0.0f, pixels, nrows, ncols, ldim, scalefactor, stridefactor, minsize, MIN(nrows, ncols));
	}

	if(!noclustering)
		ndetections = cluster_detections(rs, cs, ss, qs, ndetections);

	t = getticks() - t;

	// if the flag is set, draw each detection
	if(draw)
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=qthreshold) // check the confidence threshold
				cvCircle(frame, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 4, 8, 0); // we draw circles here since height-to-width ratio of the detected face regions is 1.0f

	// if the `verbose` flag is set, print the results to standard output
	if(verbose)
	{
		//
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=qthreshold) // check the confidence threshold
				printf("%d %d %d %f\n", (int)rs[i], (int)cs[i], (int)ss[i], qs[i]);

		//
		//printf("# %f\n", 1000.0f*t); // use '#' to ignore this line when parsing the output of the program
	}
}

void process_webcam_frames()
{
	CvCapture* capture;

	IplImage* frame;
	IplImage* framecopy;

	int stop;

	const char* windowname = "--------------------";

	//
	capture = cvCaptureFromCAM(0);

	if(!capture)
	{
		printf("* cannot initialize video capture ...\n");
		return;
	}

	// the main loop
	framecopy = 0;
	stop = 0;

	while(!stop)
	{
		// wait 5 miliseconds
		int key = cvWaitKey(5);

		// get the frame from webcam
		if(!cvGrabFrame(capture))
		{
			stop = 1;
			frame = 0;
		}
		else
			frame = cvRetrieveFrame(capture, 1);

		// we terminate the loop if the user has pressed 'q'
		if(!frame || key=='q')
			stop = 1;
		else
		{
			// we mustn't tamper with internal OpenCV buffers
			if(!framecopy)
				framecopy = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, frame->nChannels);
			cvCopy(frame, framecopy, 0);

			// webcam outputs mirrored frames (at least on my machines)
			// you can safely comment out this line if you find it unnecessary
			cvFlip(framecopy, framecopy, 1);

			// ...
			process_image(framecopy, 1);

			// ...
			cvShowImage(windowname, framecopy);
		}
	}

	// cleanup
	cvReleaseImage(&framecopy);
	cvReleaseCapture(&capture);
	cvDestroyWindow(windowname);
}

int main(int argc, char* argv[])
{
	//
	int arg;
	char input[1024], output[1024];

	//
	if(argc < 2)
	{
		printf("# please specify a valid detection cascade path\n");
		return 0;
	}
	else
	{
		int size;
		FILE* file;

		//
		file = fopen(argv[1], "rb");

		if(!file)
			return 0;

		//
		fseek(file, 0L, SEEK_END);
		size = ftell(file);
		fseek(file, 0L, SEEK_SET);

		//
		cascade = malloc(size);

		if(!cascade || size!=fread(cascade, 1, size, file))
			return 0;

		//
		fclose(file);
	}

	// set default parameters
	minsize = 128;
	maxsize = 1024;

	scalefactor = 1.1f;
	stridefactor = 0.1f;

	qthreshold = 5.0f;

	usepyr = 0;
	noclustering = 0;
	verbose = 0;

	//
	input[0] = 0;
	output[0] = 0;

	// parse command line arguments
	arg = 2;

	while(arg < argc)
	{
		//
		if(0==strcmp("-u", argv[arg]) || 0==strcmp("--usepyr", argv[arg]))
		{
			usepyr = 1;
			++arg;
		}
		else if(0==strcmp("-i", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%s", input);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-o", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%s", output);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-m", argv[arg]) || 0==strcmp("--minsize", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%d", &minsize);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-M", argv[arg]) || 0==strcmp("--maxsize", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%d", &maxsize);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-c", argv[arg]) || 0==strcmp("--scalefactor", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%f", &scalefactor);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-t", argv[arg]) || 0==strcmp("--stridefactor", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%f", &stridefactor);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-q", argv[arg]) || 0==strcmp("--qthreshold", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%f", &qthreshold);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 0;
			}
		}
		else if(0==strcmp("-n", argv[arg]) || 0==strcmp("--noclustering", argv[arg]))
		{
			noclustering = 1;
			++arg;
		}
		else if(0==strcmp("-v", argv[arg]) || 0==strcmp("--verbose", argv[arg]))
		{
			verbose = 1;
			++arg;
		}
		else
		{
			printf("# invalid command line argument '%s'\n", argv[arg]);
			return 0;
		}
	}

	if(verbose)
	{
		//
		printf("# Copyright (c) 2013, Nenad Markus\n");
		printf("# All rights reserved.\n\n");

		printf("# cascade parameters:\n");
		printf("#	tsr = %f\n", ((float*)cascade)[0]);
		printf("#	tsc = %f\n", ((float*)cascade)[1]);
		printf("#	tdepth = %d\n", ((int*)cascade)[2]);
		printf("#	ntrees = %d\n", ((int*)cascade)[3]);
		printf("# detection parameters:\n");
		printf("#	minsize = %d\n", minsize);
		printf("#	maxsize = %d\n", maxsize);
		printf("#	scalefactor = %f\n", scalefactor);
		printf("#	stridefactor = %f\n", stridefactor);
		printf("#	qthreshold = %f\n", qthreshold);
		printf("#	usepyr = %d\n", usepyr);
	}

	//
	if(0 == input[0])
		process_webcam_frames();
	else
	{
		IplImage* img;

		//
		img = cvLoadImage(input, CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("# cannot load image from '%s'\n", argv[3]);
			return 0;
		}

		process_image(img, 1);

		//
		if(0!=output[0])
			cvSaveImage(output, img, 0);

		//
		cvReleaseImage(&img);
	}

	return 0;
}
