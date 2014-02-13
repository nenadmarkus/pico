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
	$(CC) sample.c ../picort.c -lm $(OPENCV) -O3 -DQCUTOFF=3.0f -o exe

	#
	# ---------------------------

	#
	# Compiling rotation invariant face detection sample ...
	#
	$(CC) sample.c ../picorirt.c -lm $(OPENCV) -O3 -DQCUTOFF=8.0f -o exeri
