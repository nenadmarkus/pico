#
#
#

#
import sys
import os
import numpy
import cv2
import struct

#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('src')
args = parser.parse_args()

#
def write_rid_to_stdout(im):
	#
	# raw intensity data
	#

	#
	h = im.shape[0]
	w = im.shape[1]

	#
	hw = struct.pack('ii', h, w)
	pixels = struct.pack('%sB' % h*w, *im.reshape(-1))

	#
	sys.stdout.buffer.write(hw)
	sys.stdout.buffer.write(pixels)

#
for dirpath, dirnames, filenames in os.walk(args.src):
	for filename in filenames:
		if filename.endswith('.jpg') or filename.endswith('.png') or filename.endswith('.JPG') or filename.endswith('.jpeg'):
			#
			path = dirpath + '/' + filename

			#
			img = cv2.imread(path)

			if img is None:
				continue

			if len(img.shape)==3:
				img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

			#
			if max(img.shape) > 1500:
				img = cv2.resize(img, (0, 0), fx=0.5, fy=0.5)

			#
			write_rid_to_stdout(img)
			sys.stdout.buffer.write( struct.pack('i', 0) )
			#
			cv2.flip(img, 1)
			write_rid_to_stdout(img)
			sys.stdout.buffer.write( struct.pack('i', 0) )