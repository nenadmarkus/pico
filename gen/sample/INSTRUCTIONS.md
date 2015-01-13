This tutorial will guide you through the process of learning your own face detector.

## 1. Get faces (positive training samples):

* Download the GENKI dataset from <http://mplab.ucsd.edu/~nick/GENKI-R2009a.tgz>
* Extract the contents of the archive to some folder: `path/to/genki`

## 2. Get images that do not contain faces (negative training samples):

* Download a file from <https://googledrive.com/host/0Bw4IT5ZOzJj6NXlJUFh0UGZCWmc>
* Rename `0Bw4IT5ZOzJj6NXlJUFh0UGZCWmc` to `images.tar`
* Extract the contents to some folder: `path/to/nonfaces`

## 3. Start the learning process

Run the script `lrn.sh`/`lrn.bat` (depending on your OS) with appropriate parameters:

	$ lrn.* path/to/genki path/to/nonfaces

Learning should finish in ~2 days on a modern machine with 4 CPU cores.

If everything went well, you should find a classification cascade file `d` in the folder.

## 4. Try the new face detector:

* Use `picogen.c` to transform the detection cascade into highly optimized C code: `./picogen d 0.0 is_region_a_face > genki_facefinder.c`
* Modify `sample.c`: replace `#include "is_region_a_face"` with `#include "genki_facefinder.c"`
* Compile and run the program