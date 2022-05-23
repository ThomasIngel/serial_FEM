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
    cs *T;
    sed *S ;

    T = cs_load (stdin, 1) ;             /* load triplet matrix T from stdin */
    printf ("---------------------\nT triplet:\n") ;
    cs_print (T, 0) ;                    /* print T */

	S = sed_compress(T) ;                // S = sparse extr. diag. of T
	printf ("---------------------\nS sparse diag:\n") ;
	sed_print (S, 0) ;                  // print S
    
	int n = 4; 
	double b[n];
	b[0] = 10; b[1] = 30; b[2] = 65; b[3] = 119;
	
	double u[n];
	for (size_t i=0; i<n; ++i) {
		u[i] = 1;
	}

	//double tol = 1e-9;	

	cg_seriel(n, S, b, u, 1e-9);

	printf ("---------------------\n\n") ;
	printf ("TEST TEST TEST\n\n");

	printf("Vector u:\n");
    print_vec(u, n);

	cs_free (T) ;                        /* clear T */
	sed_free (S) ;                       /* clear S */
	return (0) ;
}

