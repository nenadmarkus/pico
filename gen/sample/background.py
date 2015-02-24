#
#
#

#
import sys
import os
import numpy
from PIL import Image
import struct
import argparse

#
parser = argparse.ArgumentParser()
parser.add_argument('src')
args = parser.parse_args()

#
src = args.src

#
def write_rid(im):
	#
	# raw intensity data
	#

	#
	h = im.shape[0]
	w = im.shape[1]

	#
	hw = struct.pack('ii', h, w)

	tmp = [None]*w*h
	for y in range(0, h):
		for x in range(0, w):
			tmp[y*w + x] = im[y, x]

	#
	pixels = struct.pack('%sB' % w*h, *tmp)

	#
	sys.stdout.buffer.write(hw)
	sys.stdout.buffer.write(pixels)

#
for dirpath, dirnames, filenames in os.walk(src):
	for filename in filenames:
		#
		path = dirpath + '/' + filename

		#
		try:
			im = numpy.asarray(Image.open(path).convert('L'))
		except:
			continue

		#
		write_rid(im)
		sys.stdout.buffer.write( struct.pack('i', 0) )

