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

#ifdef _WIN32
#include <opencv/cv.h>
#include <opencv/highgui.h>
#else
#include <cv.h>
#include <highgui.h>
#endif

#include "n3odet.h"

int minfacesize = 0;

int n3_find_faces(float rs[], float cs[], float ss[], float qs[], int maxndetections,
						unsigned char pixels[], int nrows, int ncols, int ldim)
{
	static char facefinder[] =
		#include "facefinder.array"
		;
	
	return n3_find_objects(rs, cs, ss, qs, maxndetections, facefinder, pixels, nrows, ncols, ldim, 1.1f, 0.1f, minfacesize, MIN(nrows, ncols), 1.5f, 1); // accurate
	///return n3_find_objects(rs, cs, ss, qs, maxndetections, facefinder, pixels, nrows, ncols, ldim, 1.2f, 0.1f, minfacesize, MIN(nrows, ncols), 0.0f, 1); // fast
}

void process_image(IplImage* frame, int draw, int print)
{
	int i;

	unsigned char* pixels;
	int nrows, ncols, ldim;

	#define MAXNDETECTIONS 2048
	int ndetections;
	float qs[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];

	static IplImage* gray = 0;

	//
	if(!gray)
		gray = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, 1);
	if(frame->nChannels == 3)
		cvCvtColor(frame, gray, CV_RGB2GRAY);
	else
		cvCopy(frame, gray, 0);

	//
	pixels = (unsigned char*)gray->imageData;
	nrows = gray->height;
	ncols = gray->width;
	ldim = gray->widthStep;

	// actually, all the smart stuff happens in the following function
	ndetections = n3_find_faces(rs, cs, ss, qs, MAXNDETECTIONS, pixels, nrows, ncols, ldim);

	// if the flag is set, draw each detection as a red circle
	if(draw)
		for(i=0; i<ndetections; ++i)
			cvCircle(frame, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 4, 8, 0);

	// if the flag is set, print the results to standard output
	if(print)
	{
		printf("%d\n", ndetections);

		for(i=0; i<ndetections; ++i)
			printf("%d %d %d\n", (int)rs[i], (int)cs[i], (int)ss[i]);
	}
}

void process_webcam_frames()
{
	CvCapture* capture;

	IplImage* frame;
	IplImage* frameCopy;
	
	int stop;

	const char* windowname = "--------------------";

	// try to initialize video capture from the default webcam
	capture = cvCaptureFromCAM(0);
	if(!capture)
	{
		printf("Cannot initialize video capture!\n");
		return;
	}

	// start the main loop in which we'll process webcam output
	frameCopy = 0;
	stop = 0;
	while(!stop)
	{
		// wait 5 miliseconds
		int key = cvWaitKey(5);

		// retrieve a pointer to the image acquired from the webcam
		if(!cvGrabFrame(capture))
			break;
		frame = cvRetrieveFrame(capture, 1);

		// we terminate the loop if we don't get any data from the webcam or the user has pressed the ESC key
		#define KEY_ESCAPE 27
		if(!frame || key==KEY_ESCAPE)
			stop = 1;
		else
		{
			// we mustn't tamper with internal OpenCV buffers and that's the reason why we're making a copy of the current frame
			if(!frameCopy)
				frameCopy = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, frame->nChannels);
			cvCopy(frame, frameCopy, 0);

			// webcam outputs mirrored frames (at least on my machines); you can safely comment out this line if you find it unnecessary
			cvFlip(frameCopy, frameCopy, 1);

			// all the smart stuff happens in the following function
			process_image(frameCopy, 1, 0);

			// display the image to the user
			cvShowImage(windowname, frameCopy);
		}
	}

	// cleanup
	cvReleaseImage(&frameCopy);
	cvReleaseCapture(&capture);
	cvDestroyWindow(windowname);
}

int main(int argc, char* argv[])
{
	IplImage* img = 0;

	if(argc==1)
	{
		/*
			The program will attempt to find faces in a video stream obtained from a default webcam attached to the computer.
			The smallest face that can be detected fits roughly in a 100x100 pixel rectangle.
		*/

		printf("Copyright (c) 2013, Nenad Markus\n");
		printf("All rights reserved.\n\n");

		minfacesize = 100;
		
		process_webcam_frames();
	}
	else if(argc==2)
	{
		/*
			The program starts by reading one integer, [MIN_FACE_SIZE], passed as a command line argument.
			The program will attempt to find faces in a video stream obtained from a default webcam attached to the computer.
			The smallest face that can be detected fits roughly in a [MIN_FACE_SIZE]x[MIN_FACE_SIZE] pixel rectangle.
		*/

		printf("Copyright (c) 2013, Nenad Markus\n");
		printf("All rights reserved.\n\n");

		sscanf(argv[1], "%d", &minfacesize);
		
		process_webcam_frames();
	}
	else if(argc==3)
	{
		/*
			The program starts by reading one integer, [MIN_FACE_SIZE], and a path to an image, [PATH], passed as command line arguments.
			The program will attempt to find faces in an image specified by [PATH].
			The smallest face that can be detected fits roughly in a [MIN_FACE_SIZE]x[MIN_FACE_SIZE] pixel rectangle.
			The program writes (to standard output) the number of detections and their positions and scales.
		*/

		sscanf(argv[1], "%d", &minfacesize);

		img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("Cannot load image!\n");
			return 1;
		}

		process_image(img, 0, 1);

		cvReleaseImage(&img);
	}
	else if(argc==4)
	{
		/*
			The program starts by reading one integer, [MIN_FACE_SIZE], path to the input image, [PATH1], and path to the output image, [PATH2], passed as command line arguments.
			The program will attempt to find faces in an image specified by [PATH1].
			The smallest face that can be detected fits roughly in a [MIN_FACE_SIZE]x[MIN_FACE_SIZE] pixel rectangle.
			The program outputs a new image to [PATH2].
			This image is just the one from [PATH1] with obtained detections drawn over it.
		*/

		sscanf(argv[1], "%d", &minfacesize);

		img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);
		if(!img)
		{
			printf("Cannot load image!\n");
			return 1;
		}

		process_image(img, 1, 1);

		cvSaveImage(argv[3], img, 0);

		cvReleaseImage(&img);
	}
	else
		return 1;

	return 0;
}
