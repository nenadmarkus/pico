#
#
#

OPENCV = -lopencv_highgui -lopencv_core -lopencv_imgproc -I/usr/include/opencv

#
#
#

facefinder:
	#
	# Compiling face detection sample ...
	#
	$(CC) sample.c ../picort.c -lm $(OPENCV) -O3 -DQCUTOFF=3.0f -o facefinder

	#
	# ---------------------------

	#
	# Compiling rotation invariant face detection sample ...
	#
	$(CC) sample.c ../picorirt.c -lm $(OPENCV) -O3 -DQCUTOFF=8.0f -o rotfacefinder

procspeed:
	#
	# Compiling processing speed tool ...
	#
	$(CC) getprocspeed.c ../picorirt.c -lm -lrt $(OPENCV) -O3 -o getprocspeed