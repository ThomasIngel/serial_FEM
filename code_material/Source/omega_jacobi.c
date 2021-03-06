#include "hpc.h"
#include "blas_level1.h"

// Serieller omega_jacobi solver

void
omega_jacobi(size_t n,
           const sed *A, const double *b, double *u, double omega, double tol,
           const double* dir, const index* dir_ind, const index n_dir) {
        // A     - stiffness matrix (sed Format!)
        // b     - righthand side
        // u     - inital guess for solution
        // omega - often 2/3
        // tol   - Toleranz (stopping criteria)
        
        // set u at dirichlet bcs to 0 because of the homogenication
        inc_dir_r(u, dir_ind, n_dir);

        double *Ax = A->x; // data of the matrix A
	
		// calculate omega * inv(D)^-1 for permormance
        double omega_inv_diag[n];
        for (index i = 0; i < n; ++i){
        	omega_inv_diag[i] = omega / Ax[i];
        }
        
        //blasl1_dcopy(Ax,diag,(index) n,1.); 
        // the diagonal of A is now in vector diag_inv	
        double r[n];
        blasl1_dcopy(b,r,(index) n,1.);         //copy b in r (r=b)

        // r = b - A*u , calculating the residuum
        sed_spmv_adapt(A,u,r,-1.0);             //Solution in vector r
        
        // residuum is zero at dirichlet nodes
        inc_dir_r(r, dir_ind, n_dir);

        // sigma = r'*r = sigma_0 , computing the scalarproduct
        double sigma_0 = blasl1_ddot(r,r,n);
        double sigma = sigma_0;

        // initializing the loop variable
        size_t k = 0;

        do {
                // increasing loop variable
                k++;

                // u_k := u_k-1 + omega * diag(A)^-1 * r
                for(index i =0;i<n;i++){
                	// compute D^-1 * r and save it in r
                    u[i] += r[i] * omega_inv_diag[i]; 
                }
                
                //blasl1_daxpy(u,r,n,omega,1.0); // u <- u + r * omega
                
                // set dirichlet nodes to 0 because of the homogenization
                inc_dir_r(u, dir_ind, n_dir);

                blasl1_dcopy(b,r,(index) n,1.);  //copy b in r (r=b)

                // r = b - A*u , calculating the residuum
                sed_spmv_adapt(A,u,r,-1.0);
                
                // residuum is 0 at dirichlet nodes
                inc_dir_r(r, dir_ind, n_dir);

                // sigma = r' * r , computing the scalarproduct with the new residuum
                sigma = blasl1_ddot(r,r,n);

        } while (sigma > tol*tol*sigma_0);
        
        // write dirichlet data in solution vector
        inc_dir_u(u, dir, dir_ind, n_dir);
}
