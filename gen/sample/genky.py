#
#
#

#
import sys
import random
import numpy
from scipy import misc
from PIL import Image
from PIL import ImageOps
import struct
import argparse
import os

#
parser = argparse.ArgumentParser()
parser.add_argument('src', help='GENKI source folder')
args = parser.parse_args()

#
src = args.src

#
plot = 0

if plot:
	import matplotlib.pyplot
	import matplotlib.image
	import matplotlib.cm

	fig = matplotlib.pyplot.figure()
	matplotlib.pyplot.show(block=False)

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
def export(im, r, c, s):
	#
	nrows = im.shape[0]
	ncols = im.shape[1]

	# crop
	r0 = max(int(r - 0.75*s), 0); r1 = min(r + 0.75*s, nrows)
	c0 = max(int(c - 0.75*s), 0); c1 = min(c + 0.75*s, ncols)

	im = im[r0:r1, c0:c1]

	nrows = im.shape[0]
	ncols = im.shape[1]

	r = r - r0
	c = c - c0

	# resize, if needed
	maxwsize = 192.0
	wsize = max(nrows, ncols)

	ratio = maxwsize/wsize

	if ratio<1.0:
		im = numpy.asarray( Image.fromarray(im).resize((int(ratio*ncols), int(ratio*nrows))) )

		r = ratio*r
		c = ratio*c
		s = ratio*s

	#
	nrands = 7;

	list = []

	for i in range(0, nrands):
		#
		stmp = s*random.uniform(0.9, 1.1)

		rtmp = r + s*random.uniform(-0.05, 0.05)
		ctmp = c + s*random.uniform(-0.05, 0.05)

		#
		if plot:
			matplotlib.pyplot.cla()

			matplotlib.pyplot.plot([ctmp-stmp/2, ctmp+stmp/2], [rtmp-stmp/2, rtmp-stmp/2], 'b', linewidth=3)
			matplotlib.pyplot.plot([ctmp+stmp/2, ctmp+stmp/2], [rtmp-stmp/2, rtmp+stmp/2], 'b', linewidth=3)
			matplotlib.pyplot.plot([ctmp+stmp/2, ctmp-stmp/2], [rtmp+stmp/2, rtmp+stmp/2], 'b', linewidth=3)
			matplotlib.pyplot.plot([ctmp-stmp/2, ctmp-stmp/2], [rtmp+stmp/2, rtmp-stmp/2], 'b', linewidth=3)

			matplotlib.pyplot.imshow(im, cmap=matplotlib.cm.Greys_r)

			matplotlib.pyplot.draw()

			response = input()

		list.append( (int(rtmp), int(ctmp), int(stmp)) )

	#
	write_rid(im)

	sys.stdout.buffer.write( struct.pack('i', nrands) )

	for i in range(0, nrands):
		sys.stdout.buffer.write( struct.pack('iii', list[i][0], list[i][1], list[i][2]) )

def mirror_and_export(im, r, c, s):
	#
	# exploit mirror symmetry of the face
	#

	# flip image
	im = numpy.asarray(ImageOps.mirror(Image.fromarray(im)))

	# flip column coordinate of the object
	c = im.shape[1] - c

	# export
	export(im, r, c, s)

# image list
imlist = open(src + '/Subsets/GENKI-SZSL/GENKI-SZSL_Images.txt', 'r').readlines()

# object sample is specified by three coordinates (row, column and size; all in pixels)
rs = [float(line.split()[1]) for line in open(src+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]
cs = [float(line.split()[0]) for line in open(src+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]
ss = [float(line.split()[2]) for line in open(src+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]

#
n = 0

for i in range(0, len(rs)):
	# construct full image path
	path = src + '/files/' + imlist[i].strip()

	r = rs[i]
	c = cs[i]
	s = ss[i]

	#
	try:
		im = Image.open(path).convert('L')
	except:
		continue

	#
	im = numpy.asarray(im)

	#
	export(im, r, c, s)

	# faces are symmetric and we exploit this here
	mirror_and_export(im, r, c, s)

