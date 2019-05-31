/*
 *  This code is released under the MIT License.
 *  Copyright (c) 2013 Nenad Markus
 */

#pragma once

/*
	
*/


// * `angle` is a number between 0 and 1 that determines the counterclockwise in-plane rotation of the cascade:
//		0.0f corresponds to 0 radians and 1.0f corresponds to 2*pi radians

int find_objects
(
	float rcsq[], int maxndetections,
	void* cascade, float angle,
	void* pixels, int nrows, int ncols, int ldim,
	float scalefactor, float stridefactor, float minsize, float maxsize
);

int cluster_detections(float rcsq[], int n);

int update_memory
(
	int* slot,
	float memory[], int counts[], int nmemslots, int maxslotsize,
	float rcsq[], int ndets, int maxndets
);