#include "hpc.h"

void get_indices_and_points(const index *Elem, const index nc, 
		const double *Coord, size_t ind[6],
		double p1[2], double p2[2], double p3[2]){
	// get indices from Elements
	for (size_t i = 0; i < 3; ++i) {
		ind[i] = Elem[i];
		ind[i+3] = nc + Elem[i+3];
	}
	
	// points from coordinates
	p1[0] = Coord[2 * ind[0]];
	p1[1] = Coord[2 * ind[0] + 1];
	p2[0] = Coord[2 * ind[1]];
	p2[1] = Coord[2 * ind[1] + 1];
	p3[0] = Coord[2 * ind[2]];
	p3[1] = Coord[2 * ind[2] + 1];
}

void calc_trans_matrix(const double p1[2], const double p2[2],
  const double p3[2], double d[2][2]){
	// loop over x and y component of the points
	for (size_t i = 0; i < 2; ++i) {
		d[0][i] = p3[i] - p2[i];
		d[1][i] = p1[i] - p3[i];
	}
}

double det_trans_matrix(const double d[2][2]) {
	// calculates the determinant of the 2 x 2 transformation matrix
	return d[0][0] * d[1][1] - d[0][1] * d[1][0];
}

void stima_local(const double p1[2], const double p2[2], const double p3[2],
  const index typ, double dx[6], double ax[9]) {
	// declare and compute transformation matrix
	double d[2][2];
	calc_trans_matrix(p1, p2, p3, d); 
	
	// compute determinant 
	double det_B = det_trans_matrix(d);
	
	// calculate scaling factor from coefficient in front of the lapalace 
	// opertor and the area of the element
	double fac = (kappa(p1,typ) + kappa(p2,typ) + kappa(p3,typ)) / 
		(6.0 * det_B);
	
	// fill the diagonal and offdiagonal entrys with the reduced computations 
	// from the lecture
	dx[0] = fac * (d[0][0] * d[0][0] + d[0][1] * d[0][1]);
	dx[1] = fac * (d[1][0] * d[1][0] + d[1][1] * d[1][1]);
	
	ax[0] = fac * (d[0][0] * d[1][0] + d[0][1] * d[1][1]);
	ax[1] = -(dx[0] + ax[0]);
	ax[2] = ax[0];
	ax[3] = -(dx[1] + ax[0]);
	ax[4] = ax[3];
	ax[5] = ax[1];
	ax[6] = 2.0 * ax[1];
	ax[7] = 2.0 * ax[3];
	ax[8] = 2.0 * ax[0];
	
	dx[2] = -(ax[1] + ax[3]);
	dx[3] = dx[0] + dx[1] + dx[2];
	dx[4] = dx[3];
	dx[5] = dx[3];
}

void diag_vec_to_sed_data(double *Tx, const size_t ind[6], const double dx[6]){
	for (size_t i = 0; i < 6; ++i) {
		Tx[ind[i]] += dx[i];
	}
}

void off_diag_vec_to_sed_data(double *Tx, const index *Ti, const size_t ind[6], 
  const double ax[9]){
  	// declare helper variables
  	size_t imin;
  	size_t imax;
  	static int ai[9] = {0, 0, 1, 1, 2, 2, 3, 3, 4}; // i indices for ax
	static int aj[9] = {3, 5, 3, 4, 4, 5, 4, 5, 5}; // j indices for ax
	
  	// loop over all off diagonal entrys
	for (size_t i = 0; i < 9; ++i) {
		// find max index
		if (ind[ai[i]] < ind[aj[i]]) {
			imin = ind[ai[i]];
			imax = ind[aj[i]];			
		} else {
			imin = ind[aj[i]];
			imax = ind[ai[i]];
		}
		
		// fill entry on right place
		for (index j = Ti[imin]; j < Ti[imin + 1]; ++j) {
			if (Ti[j] == imax) {
				Tx[j] += ax[i];
				break;
			}
		}
	}
}

index mesh_stima_global(mesh* M, sed* T) {
	// check wheter there is already data in the matrix and allocate if needed
	size_t n = T->n; // dimension of dim(T) = n x n
	index* Ti = T->i; // col pointer and row indices
	if (!(T->x)) {
		T->x = calloc(Ti[n], sizeof(double));
	}
	
	// check if allocation was sucessfull
	if (!(T->x)) {
		return(0);
	}
	
	// get some variables for readability
	size_t ne = M->nelem; // number of elements
	size_t nc = M->ncoord; // number of coordinats
	index* Elems = M->elem; // elements from the mesh
	double* Coords = M->coord; // coordinates of the nodes
	double* Tx = T->x; // pointer to matrix data
	size_t incElem = 7; // seven entrys per element
	
	
	// declare some arrays wich are needed 
	size_t ind[6]; // index array
	double p1[2]; // first point of triangle 
	double p2[2]; // second point of triangle
	double p3[2]; // third poiint of triangle 
	double dx[6]; // diagonal entries of sed matrix
	double ax[9]; // off diagonal entries of sed matrix
	
	// loop over all elements
	for (size_t k = 0; k < ne; ++k){
		// get points from the triangle and indices to construct the matrix
		get_indices_and_points(&Elems[k*incElem], nc,
			Coords, ind, 
			p1, p2, p3);
			
		// calculate local stiffnes matrix
		stima_local(p1, p2, p3, &Elems[k*incElem+7], dx, ax);
		
		// fill elements from diagonal entries
		diag_vec_to_sed_data(Tx, ind, dx);
		
		// fill elements from off diagonal entries
		off_diag_vec_to_sed_data(Tx, Ti, ind, ax);
	}
	return 1;
}
