// For Testing cg_seriell // test_cg_seriell_v2

#include "hpc.h"
#include <math.h>

void print_vec_double(const double* x, size_t n){
	for (size_t i = 0; i < n; ++i){
		printf("x[%d] = %10g \n",i,x[i]);
	}
}

void print_vec_index(const index* x, size_t n){
	for (size_t i = 0; i < n; ++i){
		printf("x[%d] = %i\n", i, x[i]);
	}
}

int
main (void) {
	// problem parameters
	size_t n = 11;
	double h = 1.0 / ((double)n - 1.0);
	index nz_max = 2 * n - 1;
	printf("h^2 = %10g\n", h*h);
	
	// initial COO matrix
	cs* T = cs_alloc(n, n, nz_max, nz_max, 0);
	
	index test;
	// fill diagonal
	for (index i = 0; i < n; ++i){
		cs_entry(T, i, i, 2 / (h * h));
	}

	// fill subdiagonal
	for (index i = 0; i < n - 1; ++i){
		cs_entry(T, i+1, i, -1/(h*h));
	}
	cs_print(T, 0);

	sed* S = sed_compress(T);
	//sed_print(S,0);
	
	// rhs and x0
	double b[n];
	double u[n];
	double u_jacobi[n];
	for (size_t i = 0; i < n; ++i){
		b[i] = 1;
		u[i] = 0;
		u_jacobi[i] = 0;
	}
	
	cg_seriell(n, S, b, u, 1e-6);
	print_vec_double(u,n);

	double omega = 2/3;
	omega_jacobi(n,S,b,u_jacobi,omega,1e-6);
	print_vec_double(u_jacobi,n);

	cs_free (T) ;                        /* clear T */
	sed_free (S) ;                       /* clear S */
	return (0) ;
}

