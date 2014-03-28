#
#
#

#
import random
import numpy
import matplotlib.pyplot
import matplotlib.image
import matplotlib.cm
from PIL import Image
from PIL import ImageOps
import struct
import argparse
import os

#
parser = argparse.ArgumentParser()
parser.add_argument('srcfolder', help='input folder, source')
parser.add_argument('dstfolder', help='output folder, destination')
args = parser.parse_args()

#
srcfolder = args.srcfolder
dstfolder = args.dstfolder

# create destination folder, if needed
if not os.path.exists(dstfolder):
   os.makedirs(dstfolder)

#
plot = 0

if plot:
	fig = matplotlib.pyplot.figure()
	matplotlib.pyplot.show(block=False)

#
def saveasrid(im, path):
	#
	# raw intensity data
	#

	#
	h = im.shape[0]
	w = im.shape[1]

	#
	f = open(path, 'wb')

	#
	data = struct.pack('ii', w, h)
	f.write(data)

	tmp = [None]*w*h
	for y in range(0, h):
		for x in range(0, w):
			tmp[y*w + x] = im[y, x]

	#
	data = struct.pack('%sB' % w*h, *tmp)
	f.write(data)

	#
	f.close()

#
def export(im, r, c, s, folder, id, list):
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
	list.write(id + '.rid\n')

	#
	nrands = 7;

	list.write('\t%d\n' % nrands)

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

			response = input('Press Enter to continue...')

		list.write('\t%f %f %f\n' % (rtmp, ctmp, stmp))

	list.write('\n')
	list.flush()

	#
	saveasrid(im, folder + '/' + id + '.rid')

def exportmirrored(im, r, c, s, folder, id, list):
	#
	# exploit mirror symmetry of the face
	#

	# flip image
	im = numpy.asarray(ImageOps.mirror(Image.fromarray(im)))

	# flip column coordinate of the object
	c = im.shape[1] - c

	# export
	export(im, r, c, s, folder, id, list)

# image list
imlist = open(srcfolder + '/Subsets/GENKI-SZSL/GENKI-SZSL_Images.txt', 'r').readlines()

# object sample is specified by three coordinates (row, column and size; all in pixels)
rs = [float(line.split()[1]) for line in open(srcfolder+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]
cs = [float(line.split()[0]) for line in open(srcfolder+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]
ss = [float(line.split()[2]) for line in open(srcfolder+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]

#
list = open(dstfolder + '/list.txt', 'w')

n = 0

for i in range(0, len(rs)):
	# image path
	path = srcfolder + '/files/' + imlist[i].strip()

	r = rs[i]
	c = cs[i]
	s = ss[i]

	#
	try:
		im = Image.open(path).convert('L')
	except:
		continue

	print(path)

	#
	im = numpy.asarray(im)

	#
	id = 'face' + str(n)
	export(im, r, c, s, dstfolder, id, list)
	n = n+1

	# faces are symmetric and we exploit this here
	id = 'face' + str(n)
	exportmirrored(im, r, c, s, dstfolder, id, list)
	n = n+1

#
list.close()
