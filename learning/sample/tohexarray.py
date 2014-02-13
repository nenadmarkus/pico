#
#
#

#
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('path')
args = parser.parse_args()

path = args.path

#
data = open(path, 'rb').read()

strdata = ''

for i in range(0, len(data)):
	strdata += hex(ord(data[i])) + ', '

print strdata