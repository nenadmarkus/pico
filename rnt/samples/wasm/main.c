#include "../../picornt.c"

#include <stdio.h>

int find_faces(
	float rcsq[],
	int maxndetections,
	unsigned char pixels[],
	int nrows,
	int ncols,
	int ldim,
	float scalefactor,
	float shiftfactor,
	int minfacesize,
	int maxfacesize
)
{
	static char facefinder[] = {
		#include "facefinder.hex"
	};

	static int slot = 0;
	static const int nmemslots = 5;
	static const int maxslotsize = 1024;
	static float memory[4*nmemslots*maxslotsize];
	static int counts[nmemslots];

	int n = find_objects(
			rcsq, maxndetections,
			facefinder,
			0.0f,
			pixels, nrows, ncols, ldim,
			scalefactor, shiftfactor,
			minfacesize, maxfacesize
	);

	n = update_memory(
		&slot,
		memory, counts, nmemslots, maxslotsize,
		rcsq, n, maxslotsize
	);

	n = cluster_detections(rcsq, n);

	return n;
}
