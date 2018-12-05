This tutorial will guide you through the process of learning your own face detector.

## 1. Get images with annotated faces:

Invoke the script `dldata.sh`:

	bash dldata.sh

This script will download the Caltech 10,000 Web Faces dataset.
This dataset contains images in which frontal faces are annotated with four facial landmark points.
These landmark points will be used to deduce the bounding box for each frontal face based on simple anthropometric relations.

The negative training samples (non-faces) are data mined from a set of image regions that do not overlap with any ground-truth bounding box.

## 2. Generate the training data

Use the script `caltechfaces.py` to convert the annotated images into the format appropriate for `picolrn`:

	python3 caltechfaces.py > trdata

The file `trdata` can now be processed with `picolrn`.

## 3. Start the training process

Start with default parameters:

	$ ./picolrn trdata d

After the learning is finished (less than a day on a modern machine), you should find a classification cascade file `d` in the folder.

## 4. Try the new face detector:

Move the detector to folder `rnt/samples/native` and execute the following command:

	$ ./pico d