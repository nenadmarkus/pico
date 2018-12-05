/*
 *  This code is released under the MIT License.
 *  Copyright (c) 2013 Nenad Markus
 */

#include <stdio.h>

// OpenCV 3.x required
// depending on your computer configuration (OpenCV install path), the following line might need modifications
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui/highgui_c.h>

//
#include "../../picornt.h"

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

float angle;

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
	float rcsq[4*MAXNDETECTIONS];

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

		ndetections = find_objects(rcsq, MAXNDETECTIONS, cascade, angle, pixels, nrows, ncols, ldim, scalefactor, stridefactor, MAX(16, minsize), MIN(128, maxsize));

		for(i=1; i<5; ++i)
		{
			cvResize(pyr[i-1], pyr[i], CV_INTER_LINEAR);

			pixels = (uint8_t*)pyr[i]->imageData;
			nrows = pyr[i]->height;
			ncols = pyr[i]->width;
			ldim = pyr[i]->widthStep;

			nd = find_objects(&rcsq[4*ndetections], MAXNDETECTIONS-ndetections, cascade, angle, pixels, nrows, ncols, ldim, scalefactor, stridefactor, MAX(64, minsize>>i), MIN(128, maxsize>>i));

			for(j=ndetections; j<ndetections+nd; ++j)
			{
				rcsq[4*j+0] = (1<<i)*rcsq[4*j+0];
				rcsq[4*j+1] = (1<<i)*rcsq[4*j+1];
				rcsq[4*j+2] = (1<<i)*rcsq[4*j+2];
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
		ndetections = find_objects(rcsq, MAXNDETECTIONS, cascade, angle, pixels, nrows, ncols, ldim, scalefactor, stridefactor, minsize, MIN(nrows, ncols));
	}

	if(!noclustering)
		ndetections = cluster_detections(rcsq, ndetections);

	t = getticks() - t;

	// if the flag is set, draw each detection
	if(draw)
		for(i=0; i<ndetections; ++i)
			if(rcsq[4*i+3]>=qthreshold) // check the confidence threshold
				cvCircle(frame, cvPoint(rcsq[4*i+1], rcsq[4*i+0]), rcsq[4*i+2]/2, CV_RGB(255, 0, 0), 2, 8, 0); // we draw circles here since height-to-width ratio of the detected face regions is 1.0f

	// if the `verbose` flag is set, print the results to standard output
	if(verbose)
	{
		//
		for(i=0; i<ndetections; ++i)
			if(rcsq[4*i+3]>=qthreshold) // check the confidence threshold
				printf("%d %d %d %f\n", (int)rcsq[4*i+0], (int)rcsq[4*i+1], (int)rcsq[4*i+2], rcsq[4*i+3]);

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
		printf("# cannot initialize video capture ...\n");
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
	if(argc<2 || 0==strcmp("-h", argv[1]) || 0==strcmp("--help", argv[1]))
	{
		printf("Usage: pico <path/to/cascade> <options>...\n");
		printf("Detect objects in images.\n");
		printf("\n");

		// command line options
		printf("Mandatory arguments to long options are mandatory for short options too.\n");
		printf("  -i,  --input=PATH          set the path to the input image\n");
		printf("                               (*.jpg, *.png, etc.)\n");
		printf("  -o,  --output=PATH         set the path to the output image\n");
		printf("                               (*.jpg, *.png, etc.)\n");
		printf("  -m,  --minsize=SIZE        sets the minimum size (in pixels) of an\n");
		printf("                               object (default is 128)\n");
		printf("  -M,  --maxsize=SIZE        sets the maximum size (in pixels) of an\n");
		printf("                               object (default is 1024)\n");
		printf("  -a,  --angle=ANGLE         cascade rotation angle:\n");
		printf("                               0.0 is 0 radians and 1.0 is 2*pi radians\n");
		printf("                               (default is 0.0)\n");
		printf("  -q,  --qthreshold=THRESH   detection quality threshold (>=0.0):\n");
		printf("                               all detections with estimated quality\n");
		printf("                               below this threshold will be discarded\n");
		printf("                               (default is 5.0)\n");
		printf("  -c,  --scalefactor=SCALE   how much to rescale the window during the\n");
		printf("                               multiscale detection process (default is 1.1)\n");
		printf("  -t,  --stridefactor=STRIDE how much to move the window between neighboring\n");
		printf("                               detections (default is 0.1, i.e., 10%%)\n");
		printf("  -u,  --usepyr              turns on the coarse image pyramid support\n");
		printf("  -n,  --noclustering        turns off detection clustering\n");
		printf("  -v,  --verbose             print details of the detection process\n");
		printf("                               to `stdout`\n");

		//
		printf("Exit status:\n");
		printf(" 0 if OK,\n");
		printf(" 1 if trouble (e.g., invalid path to input image).\n");

		//
		return 0;
	}
	else
	{
		int size;
		FILE* file;

		//
		file = fopen(argv[1], "rb");

		if(!file)
		{
			printf("# cannot read cascade from '%s'\n", argv[1]);
			return 1;
		}

		//
		fseek(file, 0L, SEEK_END);
		size = ftell(file);
		fseek(file, 0L, SEEK_SET);

		//
		cascade = malloc(size);

		if(!cascade || size!=fread(cascade, 1, size, file))
			return 1;

		//
		fclose(file);
	}

	// set default parameters
	minsize = 128;
	maxsize = 1024;

	angle = 0.0f;

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
		else if(0==strcmp("-i", argv[arg]) || 0==strcmp("--input", argv[arg]))
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
				return 1;
			}
		}
		else if(0==strcmp("-o", argv[arg]) || 0==strcmp("--output", argv[arg]))
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
				return 1;
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
				return 1;
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
				return 1;
			}
		}
		else if(0==strcmp("-a", argv[arg]) || 0==strcmp("--angle", argv[arg]))
		{
			if(arg+1 < argc)
			{
				//
				sscanf(argv[arg+1], "%f", &angle);
				arg = arg + 2;
			}
			else
			{
				printf("# missing argument after '%s'\n", argv[arg]);
				return 1;
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
				return 1;
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
				return 1;
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
				return 1;
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
			return 1;
		}
	}

	if(verbose)
	{
		//
		printf("# Copyright (c) 2013, Nenad Markus\n");
		printf("# All rights reserved.\n\n");

		printf("# cascade parameters:\n");
		printf("#	version = %d\n", ((int*)cascade)[0]);
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
			return 1;
		}

		process_image(img, 1);

		//
		if(0!=output[0])
			cvSaveImage(output, img, 0);
		else if(!verbose)
		{
			cvShowImage(input, img);
			cvWaitKey(0);
		}

		//
		cvReleaseImage(&img);
	}

	return 0;
}
