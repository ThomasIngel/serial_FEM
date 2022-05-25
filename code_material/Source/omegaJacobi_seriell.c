// omega Jacobi algorithm

#include "hpc.h"
#include "blas_level1.h"

void
cg_seriell(size_t n,
           const sed *A, const double *b, double *u, double omega, double tol) {
        // A     - stiffness matrix (sed Format!)
        // b     - righthand side
        // u     - inital guess for solution
        // omega - often 2/3
        // tol   - Toleranz (stopping criteria)

        index n = A->n; // dimension of A
        index *Ai = A->i; // index vector of the matrix A
        double *Ax = A->x; // data of the matrix A
        index nzmax = A->nzmax; // nzmax

        double diag[n];
        blasl1_dcopy(Ax,diag,(index) n,1.); // the diagonal of A is now in vector diag

        double Ax_2[nzmax];
        blasl1_dcopy(Ax,Ax_2,(index) nzmax,1.); // copy of the values of A
        for(index i =0;i<n;i++){
                Ax_2[i] -= Ax_2[i]; // replace the diagonal entries with zeros 
        }
        
        sed A_LL; // (L+L') of matrix A
        A_LL.n = n;
        A_LL.i = Ai;
        A_LL.nzmax = nzmax;
        A_LL.x = Ax_2;

        // the scalar 1-omega
        double omega_2 = 1 - omega;

        // For later calculation purposes
        double tmp[n];

        double r[n];
        blasl1_dcopy(b,r,(index) n,1.);         //copy b in r (r=b)

        // r = b - A*u , calculating the residuum
        sed_spmv_adapt(A,u,r,-1.0);             //Solution in vector r

        // sigma = r'*r = sigma_0 , computing the scalarproduct
        double sigma_0 = blasl1_ddot(r,r,n);
        double sigma = sigma_0;

        // initializing the loop variable
        size_t k = 0;

        do {
                // increasing loop variable
                k++;

                // Following computations are neede for computing:
                // A = (L + D + L') splitted into the diagonal and the upper and lower triangular matrix
                // u_k := u_k-1 + omega * D^-1 * r
                // or in a different way:
                // u_k = omega./diag(A) * (b - (L + L') * u) + (1-omega) * u
                // ./ is an element-wise division, 1-omega is saved in omega_2

                //copy b in tmp (tmp=b)
                blasl1_dcopy(b,tmp,(index) n,1.);

                // compute tmp <- b - (L + L') * u
                sed_spmv_adapt(A_LL,u,tmp,-1.0);

                // now we have left: u_k = omega * tmp/diag(A) + omega_2 * u

                // Element-wise division of tmp and diag(A)
                for(index i =0;i<n;i++){
                        tmp[i] = tmp[i]/diag[i]; 
                }

                // now we have left: u_k = omega * tmp + omega_2 * u

                blasl1_daxpy(u, tmp, n, omega, omega_2); // u <- omega_2 * tmp + omega_2 * u

                // r = b - A*u , calculating the residuum
                sed_spmv_adapt(A,u,r,-1.0);

                // sigma = r' * r , computing the scalarproduct with the new residuum
                sigma = blasl1_ddot(r,r,n);

                // Muss ich irgendwas free()?

                // Oben bei for, stimmt das oder brauche ich pointer??
               
                printf("k = %d \t norm = %10g\n", k, sqrt(sigma));
        } while (sqrt(sigma) > tol);
        
}