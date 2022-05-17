#include "hpc.h"
#include <math.h>

void vec_on_right_index(double* b, const double bx[6], const size_t ind[6]){
	for (size_t i = 0; i < 6; ++i){
		b[ind[i]] += bx[i];
	}
}

void get_midpoints(const double p1[2], const double p2[2], const double p3[2],
  double m1[2], double m2[2], double m3[2]){
	for (size_t i = 0; i < 2; ++i) {
		m1[i] = (p1[i] + p2[i]) / 2.0;
		m2[i] = (p2[i] + p3[i]) / 2.0;
		m3[i] = (p3[i] + p1[i]) / 2.0;
	} 
}

void mesh_Neumann_ti(const double p1[2], const double p2[2], const index typ, 
  const double (*fc)(double *, index), double m[3]){
	// declare some variables
	double x1[2];
	double x2[2];
	double h2 = 0.0;
	
	// calculate x1, x2 and h (not sure where this formula is coming from)
	for (size_t i = 0; i < 2; ++i){
		h2 += (p2[i] - p1[i]) * (p2[i] - p1[i]);
		x1[i] = (3.0 * p1[i] + p2[i]) / 4.0;
		x2[i] = (3.0 * p2[i] + p1[i]) / 4.0;
	}
	
	/* calculate contributions to the nodes (again not sure where this formula
	   is comming from but it is in the reference implementation)*/
	h2 = sqrt(h2) / 4.0;
	double fac = fc(x1, typ) * h2;
	m[0] = fac;
	m[2] = fac;
	
	fac = fc(x2, typ) * h2;
	m[1] = fac;
	m[2] += fac;
}

void mesh_vol_elem_ti(const double p1[2], const double p2[2], const double p3[2],
  const index typ, const double (*fc)(double *, index), double m[3]){
	// get transformation matrix and determinant
	double d[2][2];
	// for some reason here is a different used than in stima
	for (size_t i = 0; i < 2; ++i){
		d[0][i] = p1[i] - p3[i];
		d[1][i] = p2[i] - p1[i];
	}
	
	double det_B = d[0][0]*d[1][1]-d[1][0]*d[0][1];
	
	// calculate center of the triangle
	double mid[2];
	for (size_t i = 0; i < 2; ++i){
		mid[i] = (p1[i] + p2[i] + p3[i]) / 2.0;
	}
	
	// evaluate f at center, consider faktors and safe it in m
	double fac = fc(mid, typ) / 6.0 * det_B;
	for (size_t i = 0; i < 3; ++i){
		m[i] = fac;
	}
}

void mesh_vol_ti(const double p1[2], const double p2[2], const double p3[2], 
  const index typ, const double (*fc)(double *, index), double b[6]){
	// declare and calculate the midpoints
	double m1[2];
	double m2[2];
	double m3[2];
	
	get_midpoints(p1, p2, p3, m1, m2, m3);
	
	// declare vector where the contributions are safed
	double s[3];
	
	// calculate all the contributions from the virtual refined element
	mesh_vol_elem_ti(p1, m1, m3, typ, fc, s);
	b[0] = s[0];
	b[3] = s[1];
	b[5] = s[2];
	
	mesh_vol_elem_ti(p2, m2, m1, typ, fc, s);
	b[1] = s[0];
	b[4] = s[1];
	b[3] += s[3];
	
	mesh_vol_elem_ti(p3, m3, m2, typ, fc, s);
	b[2] = s[0];
	b[5] += s[1];
	b[4] += s[2];
	
	mesh_vol_elem_ti(m2, m3, m1, typ, fc, s);
	b[4] += s[0];
	b[5] += s[1];
	b[3] += s[2];
}

void mesh_RHS(const mesh* M, double* b, const double (*fV)(double *, index), 
  const double (*fN)(double *, index)){
	 // gather some variables for readability
	 size_t ne = M->nelem; // number of elements
	 size_t nc = M->ncoord; // number of coordinats
	 double* Coords = M->coord; // coordinates from mesh
	 index* Elems = M->elem; // elements from mesh
	 size_t incElem = 7; // seven entry per element
	 
	 // declare some arrays wich are needed
	 double p1[2]; // first point of triangle
	 double p2[2]; // second point of triangle
	 double p3[3]; // thrid point of triangle
	 size_t ind[6]; // indices
	 double bx[6]; // vector for contribution of the element
	 
	 // loop over all elements
	 for (size_t k = 0; k < ne; ++k) {
	 	// get indices and points from the triangle
	 	get_indices_and_points(&Elems[k*incElem], nc, Coords, ind, p1, p2, p3);
	 	
	 	// calculate volume force on the element and put in in the right place
	 	mesh_vol_ti(p1, p2, p3, Elems[k*incElem+6], fV, bx);
	 	vec_on_right_index(b, bx, ind);
	 }
	 
	 // gather some variables for readability
	 size_t nb = M->nbdry; // number of boundary conditions
	 index* Bdry = M->bdry; // bounadry conditions
	 size_t incBdry = 4; // four entries per Bdry;
	 
	 // loop over all boundarys
	 for (size_t k = 0; k < nb; ++k){
	 	if (Bdry[k*incBdry+3]){
		 	// get indices for boundary conditions
		 	for (size_t j = 0; j < 2; ++j){
		 		ind[j] = Bdry[k*incBdry + j];
		 	}
		 	// from reference implementaion no clue why
		 	ind[2] = nc + Bdry[k*incBdry + 2];
		 	
		 	// calculate contribution from neuman and save in right place
		 	// TODO: ugly form to hand over p1 and p2 write function for the points
		 	mesh_Neumann_ti(Coords + 2*ind[0], Coords + 2*ind[1], 
		 		Bdry[k*incBdry + 3], fN, bx);
		 	for (size_t j = 0; j < 3; ++j){
		 		b[ind[j]] += bx[j];
		 	}
		 }
	 }
}
