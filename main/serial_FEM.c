#include "hpc.h"
#include "blas_level1.h"
#include <errno.h>  // for errno
#include <limits.h> // for INT_MIN and INT_MAX
#include <string.h>  // for strlen
#include <stdlib.h> // for Time Bench
#include <sys/times.h>
#include <unistd.h>

// these are the functions for the boundarys, volume force etc.
double kappa( double x[2], index typ )
{
  return ( 1.0 );
}

double F_vol( double x[2], index typ )
{
  return ( 1.0 );
}

double g_Neu( double x[2], index typ )
{
        return 0.0;
  //return ( x[0] * x[1] );
}

double u_D( double x[2])
{
  return ( 2.0 );
  //return ( x[0] + x[1] );
}

double vec1_norm(const double* x, const double* y, const size_t n){
    double err = 0.0;
    double tmp;
    for (size_t i = 0; i < n; ++i){
        tmp = fabs(x[i] - y[i]);
        if (err < tmp){
            err = tmp;
                }
    }
    return err;
}

void print_vec(const double* x, size_t n){
        for (size_t i = 0; i < n; ++i){
                printf("x[%d] = %10g\n", i, x[i]);
        }
}

void print_vec2(const double* x, const double* y, size_t n){
        for (size_t i = 0; i < n; ++i){
                printf("x[%d] = %10g \t \t y[%d] = %10g\n", i, x[i], i, y[i]);
        }
}

void print_vec2_i(const double* x, const double* y, size_t n){
        for (size_t i = 0; i < n; ++i){
                printf("x[%d] = %d \t \t y[%d] = %d\n", i, x[i], i, y[i]);
        }
}

// For Benchmark
/* return real time in seconds since start of the process */
double walltime() {
    static clock_t ticks_per_second = 0;
    if (!ticks_per_second) {
        ticks_per_second = sysconf(_SC_CLK_TCK);
    }
    struct tms timebuf;
    /* times returns the number of real time ticks passed since start */
    return (double) times(&timebuf) / ticks_per_second;
}
//

int main(int argc, char** argv) {
        if (strlen(argv[1]) == 0) {
                printf("ERROR WITH REFINEMENT INPUT! ABORTING...\n");
                return 1; // empty string
        }
        char* p;
        errno = 0; // not 'int errno', because the '#include' already defined it
        long arg = strtol(argv[1], &p, 10);
        if (*p != '\0' || errno != 0) {
                printf("ERROR WITH REFINEMENT INPUT! ABORTING...\n");
                return 1; // In main(), returning non-zero means failure
        }

        if (arg < INT_MIN || arg > INT_MAX) {
                printf("ERROR WITH REFINEMENT INPUT! ABORTING...\n");
                return 1;
        }
        int norefine = arg;

        // Everything went well
        printf("Starting program with %d mesh refinement(s)\n", norefine);

        // Time
        double t0;
        double t_sm_rhs;

        t0 = walltime();

        // get mesh and stima
        mesh* H = get_refined_mesh(norefine);
        sed* A;

        A = sed_sm_build(H);
        //sed_print(A,0);

        index n = A->n;

        // get RHS
        double* b_cg = calloc(n, sizeof(double));
        double* b_jac = calloc(n, sizeof(double));
        mesh_build_rhs(H, b_cg, F_vol, g_Neu);
        mesh_build_rhs(H, b_jac, F_vol, g_Neu);


        // get number of dirichlet bcs (because of the midpoints in the mesh we
        // don't need
        index n_dir = 0;
        for (index i = 0; i < H->nfixed; ++i){
                if (H->fixed[i] >= H->ncoord){
                        break;
                }
                n_dir++;
        }

        // get dirichle boundary values and indices
        double dir[n_dir];
        get_dirich(H, u_D, dir);
        index* dir_ind = H->fixed;

        // initialize solution vector with 0
        double* u_cg = calloc(n, sizeof(double));
        double* u_jac = calloc(n, sizeof(double));

        // get dirichlet bcs for the homogenization
        for (index i = 0; i < n_dir; ++i){
                u_cg[dir_ind[i]] = dir[i];
                u_jac[dir_ind[i]] = dir[i];
        }

        // homogenitize RHS
        sed_spmv_adapt(A, u_cg, b_cg, -1.0);
        sed_spmv_adapt(A, u_jac, b_jac, -1.0);

        t_sm_rhs = walltime()-t0;
        printf("Time to build stiffness matrix & righthandside in sec: %4lf\n", t_sm_rhs);

        // Time Bench
        double t_cg;
        double t_jac;

        t0 = walltime();                                        // Time start

        double tol = 1e-6;

        // solve with cg
        cg_seriell(A, b_cg, u_cg, tol, dir, dir_ind, n_dir);
        t_cg = walltime()-t0;                                   // Time stop
        printf("Time for cg-Solver in sec: %4lf\n", t_cg);      // Time print

        // solve with \omega-jacobi
        omega_jacobi(n, A, b_jac, u_jac, 2.0 / 3.0, tol, dir, dir_ind, n_dir);
        t_jac = walltime()-t0;                                   // Time stop
        printf("Time for omega-Jacobi-Solver in sec: %4lf\n", t_jac);      // Time print

        // print solution
        //printf("cg_dir_neu \t \t jacobi_dir_neu\n");
        //print_vec2(u_cg, u_jac, n);

        // free allocated memory
        free(b_cg);
        free(b_jac);
        free(u_cg);
        free(u_jac);

        return 0;
}
