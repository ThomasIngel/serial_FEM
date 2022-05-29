// cg_seriell Version2

#include "hpc.h"
#include "blas_level1.h"

void inc_dir_res(double* r, const index* dir_ind, const size_t n_dir){
	for (size_t i = 0; i < n_dir; ++i){
		r[dir_ind[i]] = 0;
	}
}

void inc_dir_u(double* u, const double* dir, const index* dir_ind, 
  const size_t n_dir){
	for (size_t i = 0; i < n_dir; ++i){
		u[dir_ind[i]] = dir[i];
	}
}

void print_vec(const double* x, const size_t n){
	for (size_t i = 0; i < n; ++i){
		printf("x[%d] = %10g\n", i, x[i]);
	}
}

void
cg_seriell(const size_t n,
           const sed *A, const double *b, double *u, const double tol,
           const double* dir, const index* dir_ind, const size_t n_dir) {
        // A   - stiffness matrix (sed Format!)
        // b   - righthand side
        // u   - inital guess for solution
        // tol - Toleranz (stopping criteria)

		inc_dir_u(u, dir, dir_ind, n_dir);
        double r[n];
        blasl1_dcopy(b,r,(index) n,1.);         //kopiert b in r (also r=b)

        // r = b - A*u    r = r-A*u
        sed_spmv_adapt(A,u,r,-1.0);             //Ergebnis steht in r
        inc_dir_res(r, dir_ind, n_dir);
        
        //sed_gaxpy(A, u, zeros);
        //blasl1_daxpy(r, zeros, (index) n, -1.0, 1.0);

        // sigma = r'*r (Skalarprodukt)
        double sigma_0 = blasl1_ddot(r,r,n);
        double sigma = sigma_0;

        // d = r
        double d[n];
        blasl1_dcopy(r,d,(index) n,1.0);

        size_t k = 0;
        do {
                k++;
                // ad = A*d (damit nur 1x Matrixprodukt)
                double *ad = calloc(n, sizeof(double));         //statt double ad[n]    ad mit 0en initiieren calloc
                if (ad == NULL) {
                        printf("Error! memory not allocated.");
                        exit(0);
                }
                sed_spmv_adapt(A,d,ad,1.);
                
                //sed_gaxpy(A, d, ad);

                // alpha = sigma/(d*ad)
                double dad = blasl1_ddot(d,ad,n);
                double alpha = sigma / dad;

                // Update: u = u + alpha*d
                blasl1_daxpy(u, d, (index) n, alpha, 1.0);
                inc_dir_u(u, dir, dir_ind, n_dir);

                // r = r - alpha*ad
                blasl1_daxpy(r, ad, (index) n, -alpha, 1.0);
                inc_dir_res(r, dir_ind, n_dir);

                // sigma_neu = r' * r
                double sigma_neu = blasl1_ddot(r,r,n);

                // d = (sigma_neu/sigma)*d + r
                blasl1_daxpy(d, r, (index) n, 1.0, sigma_neu / sigma);

                // sigma = sigma_neu
                sigma = sigma_neu;
                free(ad);
              printf("k = %d \t norm = %10g\n", k, sqrt(sigma));
        } while (sqrt(sigma) > tol);
        //free(zeros);
}
