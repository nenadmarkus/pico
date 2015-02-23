#
#
#

#
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
def print_as_rid(im):
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
	sys.stdout.write(hw)
	sys.stdout.write(pixels)

#
for dirpath, dirnames, filenames in os.walk(srcfolder):
	for filename in filenames:
		#
		path = dirpath + '/' + filename

		#
		try:
			im = numpy.asarray(Image.open(path).convert('L'))
		except:
			continue

		#
		print_as_rid(im)
		sys.stdout.write( struct.pack('i', 0) )
