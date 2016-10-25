This tutorial will guide you through the process of learning your own face detector.

## 1. Get images with annotated faces:

* Download the GENKI dataset from <http://mplab.ucsd.edu/~nick/GENKI-R2009a.tgz> (or from <https://docs.google.com/uc?export=download&confirm=cx-l&id=0Bw4IT5ZOzJj6bm5scnhQWlltaUE>, if the original link is not working)
* Extract the contents of the archive to some folder: `path/to/genki`

## 2. Get images that do not contain faces (background):

* Download `background.tar` from <https://drive.google.com/file/d/0Bw4IT5ZOzJj6Y2ZjRzRMMjNkLTQ/view?usp=sharing>
* Extract the contents to some folder: `path/to/background`

The negative training samples are data mined from these images during the training process.

## 3. Generate the training data

* Prepare positive training samples: `python genki.py path/to/genki > trdata`
* Prepare background images: `python background.py path/to/background >> trdata`

(Important: both scripts are written in Python 3!)

The file `trdata` can now be processed with `picolrn`.

## 4. Start the training process

Start with default parameters:

	$ ./picolrn trdata d

After the learning is finished (~2 days on a modern machine with 4 CPU cores), you should find a classification cascade file `d` in the folder.

## 5. Try the new face detector:

Move the detector to folder `rnt/sample/` and execute the following command:

	$ ./pico d
