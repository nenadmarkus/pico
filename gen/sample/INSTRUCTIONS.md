This tutorial will guide you through the process of learning your own face detector.

## 1. Get images with annotated faces:

* Download the GENKI dataset from <http://mplab.ucsd.edu/~nick/GENKI-R2009a.tgz>
* Extract the contents of the archive to some folder: `path/to/genki`

## 2. Get images that do not contain faces (background):

* Download a file from <https://googledrive.com/host/0Bw4IT5ZOzJj6NXlJUFh0UGZCWmc>
* Rename `0Bw4IT5ZOzJj6NXlJUFh0UGZCWmc` to `background.tar`
* Extract the contents to some folder: `path/to/background`

The negative training samples are data mined from these images during the training process.

## 3. Generate the training data

* Prepare positive training samples: `python genki.py path/to/genki > trdata`
* Prepare background images: `python background.py path/to/background >> trdata`

The file `trdata` can now be processed with `picolrn`.

## 4. Start the training process

Start with default parameters:

	$ ./picolrn trdata d

After the learning is finished (~2 days on a modern machine with 4 CPU cores), you should find a classification cascade file `d` in the folder.

## 5. Try the new face detector:

* Use `picogen.c` to transform the detection cascade into highly optimized C code: `./picogen d 0.0 run_facefinder > facefinder_genki.c`
* Modify `sample.c`: replace `#include "facefinder.c"` with `#include "facefinder_genki.c"`
* Compile and run the program