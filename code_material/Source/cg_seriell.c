#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>

#include "blas_level1.c"
#include "sed_spmv.c"

void
cg_seriell(size_t n,
	   const sed *A, const double *b, double *u, double tol) {
	// A   - stiffness matrix (sed Format!)
	// b   - righthand side
	// u   - inital guess for solution
	// tol - Toleranz (stopping criteria)
	
	double r[n];
	blasl1_icopy(b,r,n,1.);		//kopiert b in r (also r=b)
	
	// r = b - A*u    r = r-A*b
	sed_spmv(A,u,r,-1.,1.);		//Ergebnis steht in r

	// sigma = r'*r (Skalarprodukt)
	double sigma_0 = blasl1_ddot(r,r,n);
	double sigma = sigma_0;

	// d = r
	double d[n];
	blasl1_icopy(r,d,n,1.);	

	do {
		// ad = A*d (damit nur 1x Matrixprodukt)
		double ad[n];	
		sed_spmv(A,d,ad,1.,0.);

		// alpha = sigma/(d*ad)
		double dad = blasl1_ddot(d,ad,n);
		double alpha = sigma/dad;
			
		// Update: u = u + alpha*d
		blasl1_daxpy(u,d,(index) n,alpha,1.);

		// r = r - alpha*ad
		blasl1_daxpy(r,ad,(index) n,-alpha,1.);		

		// sigma_neu = r' * r
		double sigma_neu = blasl1_ddot(r,r,n);

		// d = (sigma_neu/sigma)*d + r
		blasl1_daxpy(d,r,(index) n,1.,sigma_neu/sigma);
		
		// sigma = sigma_neu
		sigma = sigma_neu;

	} while (sqrt(sigma/sigma_0) > tol);
		
}
