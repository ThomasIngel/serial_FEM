#include "hpc.h"

void get_dirich(const mesh* H, double (*f_dir)(double *), double* dir){
	// gather some variables for readability
	index nfixed = H->nfixed;
	index* fixed = H->fixed;
	double* coord = H->coord;
	index ncoord = H->ncoord;
	
	// declare some helper variable
	double x[2];
	for (index i = 0; i < nfixed; ++i){
		// because we have the midpoints we brake after ncoord (fixed is sorted
		// and midoints are located after ncoords); we dont need this for the
		// local mesh struct
		if (fixed[i] >= ncoord){
			break;
		}
		// get coordinates of points
		x[0] = coord[2 * fixed[i]];
		x[1] = coord[2 * fixed[i] + 1];
		
		// safe dirichlet bc in right place
		dir[i] = f_dir(x);
	}
}
