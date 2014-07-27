#!/bin/bash

#
# prepare face samples (from the GENKI dataset)
#

mkdir -p faces

python preparefacesamplesfromgenki.py $1 faces

#
# prepare non-face samples (background)
#

mkdir -p nonfaces

python preparebackground.py $2 nonfaces

#
# start the learning process
#

# create an object detector
./picolrn 1 1 d > log.txt

# append stages
./picolrn d faces nonfaces 1 1e-6 6 0.980 0.5 1 d >> log.txt
./picolrn d faces nonfaces 1 1e-6 6 0.985 0.5 1 d >> log.txt
./picolrn d faces nonfaces 1 1e-6 6 0.990 0.5 2 d >> log.txt
./picolrn d faces nonfaces 1 1e-6 6 0.995 0.5 3 d >> log.txt
./picolrn d faces nonfaces 6 1e-6 6 0.997 0.5 10 d >> log.txt
./picolrn d faces nonfaces 10 1e-6 6 0.999 0.5 20 d >> log.txt

#
# transform detector to a hex array
#

python tohexarray.py d > face-detector-from-genki-dataset.ea
