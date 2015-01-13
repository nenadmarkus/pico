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
	$(CC) sample.c ../picornt.c -lm -lrt $(OPENCV) -O3 -DQCUTOFF=3.0 -o pico
