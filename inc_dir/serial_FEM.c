#include "hpc.h"

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
  //return ( x[0] * x[1] );
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
   	
   	double* b = calloc(n, sizeof(double));
   	mesh_build_rhs(H, b, F_vol, g_Neu);
   	
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
   	
   	double* u = calloc(n, sizeof(double));
   	print_vec(b,n);
   	cg_seriell(A, b, u, 1e-6, dir, dir_ind, n_dir);
   	print_vec(u,n);
   	
   	free(b);
   	free(u);
   	
   	
   	
   	
   	
   	
   	
   	
   	 
}
