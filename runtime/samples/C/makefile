#
#
#

OPENCV = -lopencv_highgui -lopencv_core -lopencv_imgproc -I/usr/include/opencv

#
#
#

output:
	#
	# Compiling face detection sample ...
	#
	$(CC) sample.c -lm $(OPENCV) -O3 -DQCUTOFF=3.0 -o exe

utils:
	#
	# Compiling processing speed tool ...
	#
	$(CC) procspeed.c -lm -lrt $(OPENCV) -O3 -o procspeed
