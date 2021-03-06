#include "hpc.h"
#include "blas_level1.h"

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
  return ( 0.0 );
  return ( x[0] + x[1] );
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

int main() {
	// inital mesh
    mesh* H = get_refined_mesh(1);
    sed* A;
    
   	A = sed_sm_build(H);
   	sed_print(A,0);
   	
   	index n = A->n;
   	
   	double* b_cg = calloc(n, sizeof(double));
   	double* b_jac = calloc(n, sizeof(double));
   	
   	mesh_build_rhs(H, b_cg, F_vol, g_Neu);
   	mesh_build_rhs(H, b_jac, F_vol, g_Neu);
   	
   	printf("b wo dir = \n");
   	print_vec(b_cg, n);
   	
   	// get dirichlet bcs
   	index n_dir = 0;
   	for (index i = 0; i < H->nfixed; ++i){
   		if (H->fixed[i] >= H->ncoord){
   			break;
   		}
   		n_dir++;
   	}
   	
   	double dir[n_dir];
   	get_dirich(H, u_D, dir);
   	
   	index* dir_ind = H->fixed;
   	
   	double* u_cg = calloc(n, sizeof(double));
   	double* u_jac = calloc(n, sizeof(double));
   	
   	for (index i = 0; i < n_dir; ++i){
   		u_cg[dir_ind[i]] = dir[i];
   		u_jac[dir_ind[i]] = dir[i];
   	}
   	
   	print_vec(u_cg, n);
   	sed_spmv_adapt(A, u_cg, b_cg, -1.0);
   	sed_spmv_adapt(A, u_jac, b_jac, -1.0);
   	
   	printf("b_cg = \n");
   	print_vec(b_cg,n);
   	cg_seriell(A, b_cg, u_cg, 1e-6, dir, dir_ind, n_dir);
   	print_vec(u_cg,n);
   	printf("-----------------\n");
   	
   	
   	omega_jacobi(n, A, b_jac, u_jac, 2.0 / 3.0, 1e-6, dir, dir_ind, n_dir);
   	print_vec(u_jac, n);
   	
   	free(b_cg);
   	free(b_jac);
   	free(u_cg);
   	free(u_jac);
}
