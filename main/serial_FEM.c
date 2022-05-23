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
	//return 0.0;
  return ( x[0] * x[1] );
}

double u_D( double x[2])
{
//  return ( 0.0 );
  return ( x[0] * x[1] );
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
		printf("x[%d] = %10g\n", i+1, x[i]);
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
    sed* A = sed_nz_pattern(H);
    sed* B = sed_nz_pattern(H);
    
    // construct my stiffnes matrix and the reference one
    mesh_stima_global(H, A);
    sed_buildS(H, B);
    //sed_print(A, 0);
    //sed_print(B, 0);
    
    size_t n = A->n;
    // my RHS and the reference one
    double* b1 = calloc(n, sizeof(double));
    double* b2 = calloc(n, sizeof(double));
    double* u = calloc(n, sizeof(double));
    
    mesh_RHS(H, b1, F_vol, g_Neu);
    mesh_buildRhs(H, b2, F_vol, g_Neu);

    // check matrix and RHS
    double err_A = vec1_norm(A->x, B->x, A->nzmax);
    printf("Error of matrix is %10g\n", err_A);

    double err_b = vec1_norm(b1, b2, n);
    printf("Error of RHS is %10g\n", err_b);
    
    double* y1 = calloc(n, sizeof(double));
    double* y2 = calloc(n, sizeof(double));
	sed_spmv_adapt(A, b1, y1, 1.0);	
	sed_gaxpy(B, b2, y2);
	
	double err_spmv = vec1_norm(y1, y2, n);
	printf("err spmv = %10g\n", err_spmv);
    // TODO: Solve LSE
    cg_seriel(n, B, b2, u, 1e-6);
    
    free(b1);
    free(b2);
    free(u);
    free(y1);
    free(y2);
}
