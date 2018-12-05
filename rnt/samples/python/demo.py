import subprocess
import numpy
import os
import ctypes
import cv2
import time

#
#
#

os.system('cc ../../picornt.c -O3 -fPIC -shared -o picornt.lib.so')
pico = ctypes.cdll.LoadLibrary('./picornt.lib.so')
os.system('rm picornt.lib.so')

#
#
#

bytes = open('../../cascades/facefinder', 'rb').read()
cascade = numpy.frombuffer(bytes, dtype=numpy.uint8)

def process_frame(gray):
	#
	maxndets = 2048
	dets = numpy.zeros(4*maxndets, dtype=numpy.float32)

	ndets = pico.find_objects(
		ctypes.c_void_p(dets.ctypes.data), ctypes.c_int(maxndets),
		ctypes.c_void_p(cascade.ctypes.data), ctypes.c_float(0.0),
		ctypes.c_void_p(gray.ctypes.data), ctypes.c_int(gray.shape[0]), ctypes.c_int(gray.shape[1]), ctypes.c_int(gray.shape[1]),
		ctypes.c_float(1.1), ctypes.c_float(0.1), ctypes.c_float(100), ctypes.c_float(1000)
	)

	ndets = pico.cluster_detections(
		ctypes.c_void_p(dets.ctypes.data), ctypes.c_int(ndets)
	)

	return list(dets.reshape(-1, 4))[0:ndets]

#
#
#

cap = cv2.VideoCapture(0)

while(True):
	#
	ret, frm = cap.read()
	#frm = cv2.resize(frm, (0,0), fx=0.5, fy=0.5)
	#
	gray = numpy.ascontiguousarray(frm[:, :, 1].reshape((frm.shape[0], frm.shape[1])))
	#
	dets = process_frame(gray) # gray needs to be numpy.uint8 array
	for det in dets:
		if det[3] >= 5.0:
			cv2.circle(frm, (int(det[1]), int(det[0])), int(det[2]/2.0), (0, 0, 255), 4)
	#
	cv2.imshow('...', frm)
	#
	if cv2.waitKey(5) & 0xFF == ord('q'):
		break

cap.release()
cv2.destroyAllWindows()