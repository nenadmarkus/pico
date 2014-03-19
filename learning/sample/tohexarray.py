#
#
#

#
import argparse

#
parser = argparse.ArgumentParser()
parser.add_argument('path')
args = parser.parse_args()

path = args.path

#
data = open(path, 'rb').read()

if len(data)!=0:
	#
	#
	#

	#
	strdata = ''

	for i in range(0, len(data)-1):
		strdata += hex(data[i]) + ', '
	strdata += hex(data[len(data)-1])

	print(strdata)