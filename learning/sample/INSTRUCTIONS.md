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

If everything went well, you should find a file `face-detector-from-genki-dataset.ea` in the folder.

## 4. Try the new face detector:

* Open `sample.c` with a text editor and replace the line `facefinder.ea` with `face-detector-from-genki-dataset.ea`
* Compile and run the program
