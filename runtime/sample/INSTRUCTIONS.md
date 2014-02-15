## Compiling the program

To compile and run the programs, you will need:

* a C compiler
* OpenCV library, <http://opencv.org/>

The programs can be compiled on *nix machines by using the provided makefile (some small modifications are usually required).

The programs use OpenCV library to obtain data from webcams and to load/save images in common image formats (JPG, PNG, ...).
Thus, if you wish to redistribute them or their modifications in binary form, they have to reproduce the OpenCV license (read more at the official web page).

## Invoking the programs

There are four ways of invoking the face detection sample program:

1. Run the program without any arguments (this is equivalent to double-clicking the executable file). In this case, the program will attempt to find faces in a video stream obtained from a default webcam attached to the computer. The smallest face that can be detected fits roughly in a 100x100 pixel rectangle.

		$ ./exe

2. Run the program by passing one integer, `MINFACESIZE`, as a command line argument. The program will attempt to find faces in a video stream obtained from a default webcam attached to the computer. The smallest face that can be detected fits roughly in a MINFACESIZExMINFACESIZE pixel rectangle.

		$ ./exe 50

3. Run the program by passing one integer, `MINFACESIZE`, and path to the input image, `PATH`, as command line arguments. The program will attempt to find faces in an image specified by `PATH`. The smallest face that can be detected fits roughly in a `MINFACESIZE x MINFACESIZE` pixel rectangle. The program writes (to standard output) the number of detections and their positions and scales.

		$ ./exe 50 /some-folder/image.jpg

4. Run the program by passing one integer, `MINFACESIZE`, path to the input image, `PATH1`, and path to the output image, `PATH2`, as command line arguments. The program will attempt to find faces in the image specified by `PATH1`. The smallest face that can be detected fits roughly in a `MINFACESIZE x MINFACESIZE` pixel rectangle. The program outputs a new image to `PATH2`. This image is just the one from `PATH1` with obtained detections drawn over it.

		$ ./exe 50 /some-folder/input-image.jpg /other-folder/output-image.png
