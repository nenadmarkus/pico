# Pixel Intensity Comparison-based Object detection (pico)

This repository contains the code for real-time face detection.
Check out a demo video at <http://www.youtube.com/watch?v=1lXfm-PZz0Q> to get the better idea.

The **pico** framework is a modifcation of the standard Viola-Jones method.
The basic idea is to scan the image with a cascade of binary classifers at all reasonable positions and scales.
An image region is classifed as an object of interest if it successfully passes all the members of the cascade.
Each binary classifier consists of an ensemble of decision trees with pixel intensity comparisons as binary tests in their internal nodes.
This enables the detector to process image regions at very high speed.
The details are given in <http://arxiv.org/abs/1305.4537>.

Some highlights of **pico** are:

* high processing speed;
* there is no need for image preprocessing prior to detection;
* there is no need for the computation of integral images, image pyramid, HOG pyramid or any other similar data structure;
* all binary tests in internal nodes of the trees are based on the same feature type (not the case in the V-J framework);
* the method can easily be modified for fast detection of in-plane rotated objects.

It can be said that the main limitation of **pico** is also its simplicity:
**pico** should be avoided when there is a large variation in appearance of the object class.
This means, for example, that **pico** should not be used for detecting pedestrians.
Large, modern convolutional neural networks are suitable for such cases.

However, **pico** can be used for simple object classes (e.g., faces or templates) when real-time performance is desired.

## Detecting objects in images and videos

The folder **rnt/** contains all the needed resources to perform object detection in images and video streams with pre-trained classification cascades.
Specifically, sample applications that perform face detection can be found in the folder **rnt/samples/**.

Note that the library also enables the detection of rotated objects without the need of image resampling or classification cascade retraining.
This is achieved by rotating the binary tests in internal tree nodes, as described in the paper.

### Embedding the runtime within your application

To use the runtime in your own application, you have to:

* Include a prototype of `find_objects(...)` in your code (for example, by adding `#include picornt.h`)
* Compile `picornt.c` with your code
* Load/include a classification cascade generated with `picolrn` (e.g., `rnt/cascades/facefinder`)
* Invoke `find_objects(...)` with appropriate parameters

Notice that there are no specific library dependencies, i.e., the code can be compiled out-of-the-box with a standard C compiler.

To get a feel for how the library works, we recommend that you look at `rnt/samples/native/sample.c` as is was specifically written to be used as documentation.

## Learning custom object detectors

The program `picolrn.c` (available in the folder **gen/**) enables you to learn your own (custom) object detectors.
The training data has to be provided in a specific format.
The details are printed to the standard output when `picolrn` is invoked without parameters.
It is often convenient to pipe this information into a text file:

    $ ./picolrn > howto.txt

A tutorial that guides you through the process of learning a face detector can be found in the folder **gen/sample/**.

## Citation

If you use the provided code/binaries for your work, please cite the following paper:
> N. Markus, M. Frljak, I. S. Pandzic, J. Ahlberg and R. Forchheimer, "Object Detection with Pixel Intensity Comparisons Organized in Decision Trees", <http://arxiv.org/abs/1305.4537>

## Contact

For any additional information contact me at <nenad.markus@fer.hr>.

Copyright (c) 2013, Nenad Markus.
All rights reserved.