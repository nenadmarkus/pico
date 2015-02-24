#!/bin/bash

#
# prepare face samples (from the GENKI dataset)
#

python genky.py $1 > trdata

#
# prepare non-face samples (background)
#

python background.py $2 >> trdata

#
# start the learning process
#

# create an object detector
./picolrn 1 1 6 d > log.txt

# append stages
./picolrn d trdata 1 1e-6 0.980 0.5 1 d >> log.txt
./picolrn d trdata 1 1e-6 0.985 0.5 1 d >> log.txt
./picolrn d trdata 1 1e-6 0.990 0.5 2 d >> log.txt
./picolrn d trdata 1 1e-6 0.995 0.5 3 d >> log.txt
./picolrn d trdata 6 1e-6 0.997 0.5 10 d >> log.txt
./picolrn d trdata 10 1e-6 0.999 0.5 20 d >> log.txt
