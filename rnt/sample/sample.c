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

#include "facefinder.c"

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
	
*/

int minsize = 100;

void process_image(IplImage* frame, int draw, int print)
{
	int i;
	float t;

	uint8_t* pixels;
	int nrows, ncols, ldim;

	#define MAXNDETECTIONS 2048
	int ndetections;
	float qs[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];

	static IplImage* gray = 0;

	/*
		IMPORTANT: these parameters are highly specific for each detection cascade and should be determined experimentally
	*/

	int (*run_detection_cascade)(float*, int, int, int, uint8_t*, int, int, int) = run_facefinder;

	float qthreshold = 5.0f; // detection quality threshold, you can vary TPR and FPR with this value
	float scalefactor = 1.1f; // how much to rescale the window during the multiscale detection process
	float stridefactor = 0.1f; // how much to move the window between neighboring detections

	/*
		...
	*/

	// grayscale image
	if(!gray)
		gray = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, 1);
	if(frame->nChannels == 3)
		cvCvtColor(frame, gray, CV_RGB2GRAY);
	else
		cvCopy(frame, gray, 0);

	// get relevant image data
	pixels = (uint8_t*)gray->imageData;
	nrows = gray->height;
	ncols = gray->width;
	ldim = gray->widthStep;

	// perform detection with the pico library
	t = getticks();
	ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, run_detection_cascade, pixels, nrows, ncols, ldim, scalefactor, stridefactor, minsize, MIN(nrows, ncols));
	ndetections = cluster_detections(rs, cs, ss, qs, ndetections);
	t = getticks() - t;

	// if the flag is set, draw each detection
	if(draw)
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=qthreshold) // check the confidence threshold
				cvCircle(frame, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 4, 8, 0); // we draw circles here since height-to-width ratio of the detected face regions is 1.0f

	// if the flag is set, print the results to standard output
	if(print)
	{
		//
		for(i=0; i<ndetections; ++i)
			if(qs[i]>=qthreshold) // check the confidence threshold
				printf("%d %d %d %f\n", (int)rs[i], (int)cs[i], (int)ss[i], qs[i]);

		//
		printf("# %f\n", 1000.0f*t); // use '#' to ignore this line when parsing the output of the program
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
			process_image(framecopy, 1, 0);

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
	IplImage* img = 0;

	if(argc==1)
	{
		printf("Copyright (c) 2013, Nenad Markus\n");
		printf("All rights reserved.\n\n");

		process_webcam_frames();
	}
	else if(argc==2)
	{
		printf("Copyright (c) 2013, Nenad Markus\n");
		printf("All rights reserved.\n\n");

		sscanf(argv[1], "%d", &minsize);

		process_webcam_frames();
	}
	else if(argc==3)
	{
		sscanf(argv[1], "%d", &minsize);

		img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("* cannot load image!\n");
			return 1;
		}

		process_image(img, 0, 1);

		cvReleaseImage(&img);
	}
	else if(argc==4)
	{
		sscanf(argv[1], "%d", &minsize);

		img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("* cannot load image!\n");
			return 1;
		}

		process_image(img, 1, 0);

		//
		cvSaveImage(argv[3], img, 0);

		//
		cvReleaseImage(&img);
	}
	else
		return 1;

	return 0;
}
