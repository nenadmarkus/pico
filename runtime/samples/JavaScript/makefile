#
#
#

#
CC=/path/to/emcc

#
output:
	$(CC) -O2 -s WARN_ON_UNDEFINED_SYMBOLS=1 link.c -o sample.html -s EXPORTED_FUNCTIONS="['_find_faces', '_malloc', '_free']" --post-js post.js
