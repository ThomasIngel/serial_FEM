// For Testing cg_seriell // test_cg_seriell_v2

#include "hpc.h"
#include <math.h>

void print_vec(double* x, size_t n){
	for (size_t i = 0; i < n; ++i){
		printf("x[%d] = %10g \n",i,x[i]);
	}
}

int
main (void) {
	// problem parameters
	size_t n = 11;
	double h = (1 / (n - 1)) * (1 / (n - 1));

	// initial COO matrix
	cs* T;
	double* data = T->x;
	index* ind_row = T->ind;
	index* ind_col = T->p;
	T->nzmax = 2 * n - 1;
	T->nz = 2 * n - 1;
	T->m = n;
	T->n = n;

	// first entrys
	data[0] = 2 * h;
	ind_row[0] = 0;
	ind_col[0] = 0;
	// last entrys
	data[nz] = 2 * h;
	ind_row[nz] = nz;
	ind_col[nz] = nz;

	for (size_t i = 1; i < n - 1; ++i){
		// diagonal
		data[i] = 2 * h;
		ind_row[i] = i;
		ind_col[i] = i;

		// subdiagonal
		data[n+i] = -h;
		ind_row[n+i] = i + 1;
		ind_col[n+i] = i;
	}

	sed* S = sed_compress(T);
	sed_print(S,0);

	cs_free (T) ;                        /* clear T */
	sed_free (S) ;                       /* clear S */
	return (0) ;
}

