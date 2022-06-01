// cg_seriell Version 27.05.22

#include "hpc.h"
#include "blas_level1.h"

void inc_dir_u(double* u, const double* dir, const index* dir_ind, 
		const index n_dir){
	for (index i = 0; i < n_dir; ++i){
		u[dir_ind[i]] = dir[i];
		printf("dir_ind[%d] = %d\n", i, dir_ind[i]);
	}		
}

void inc_dir_r(double* r, const index* dir_ind, const index n_dir){
	for (index i = 0; i < n_dir; ++i){
		r[dir_ind[i]] = 0;
	}
}

void
cg_seriell(const sed *A, const double *b, double *u, double tol,
		const double* dir, const index* dir_ind, const index n_dir) {
        // A   - stiffness matrix (sed Format!)
        // b   - righthand side
        // u   - inital guess for solution
        // tol - Toleranz (stopping criteria)
        
        // incoporate dirichlet bcs into u
        inc_dir_u(u, dir, dir_ind, n_dir);

        index n = A->n ;                         //Matrix Dim
        print_vec(u, n);

        double r[n];
        blasl1_dcopy(b,r,n,1.0);         //kopiert b in r (also r=b)


        // r = b - A*u    r = r-A*u
        /*TestBeta: sed_spmv_adapt(A,u,r,-1.0,1);             //Ergebnis steht in r */
        sed_spmv_adapt(A,u,r,-1.0);
        inc_dir_r(r, dir_ind, n_dir);

        // sigma = r'*r (Skalarprodukt)
        double sigma_0 = blasl1_ddot(r,r, (size_t) n);
        double sigma = sigma_0;

        // d = r
        double d[n];
        blasl1_dcopy(r,d,n,1.0);

        // Speicher allokieren für ad
        double *ad = calloc(n, sizeof(double));         // ad mit 0en initiieren calloc
        if (ad == NULL) {
                printf("Error! memory not allocated.");
                exit(0);
        }

        size_t k = 0;
         do {
                k++;

                // ad = A*d (damit nur 1x Matrixprodukt)
                if (k>0) {                              // ad mit 0en initiieren
                    for (index i=0; i<n; i++) {
                        ad[i] = 0;
                    }
                }
                // ad = A*d (damit nur 1x Matrixprodukt)
                //TestBeta: sed_spmv_adapt(A,d,ad,1.0,0);
                sed_spmv_adapt(A,d,ad,1.0);

                // alpha = sigma/(d*ad)
                double dad = blasl1_ddot(d,ad, (size_t) n);
                double alpha = sigma / dad;

                // Update: u = u + alpha*d
                blasl1_daxpy(u, d, n, alpha, 1.0);

                // r = r - alpha*ad
                blasl1_daxpy(r, ad, n, -alpha, 1.0);
                inc_dir_r(r, dir_ind, n_dir);

                // sigma_neu = r' * r
                double sigma_neu = blasl1_ddot(r,r, (size_t) n);

                // d = (sigma_neu/sigma)*d + r
                blasl1_daxpy(d, r, n, 1.0, sigma_neu / sigma);

                // sigma = sigma_neu
                sigma = sigma_neu;
                /* free(ad); */

                printf("k = %d \t norm = %10g\n", k, sqrt(sigma));

        } while (sqrt(sigma) > tol);
        free(ad);
}
