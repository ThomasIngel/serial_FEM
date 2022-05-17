#include "hpc.h"

// these are the functions for the boundarys, volume force etc.
double kappa( double x[2], index typ )
{
  return ( 1.0 );
}

double F_vol( double x[2], index typ )
{
  return ( 0.0 );
}

double g_Neu( double x[2], index typ )
{
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
        tmp = abs(x[i] - y[i]);
        if (err < tmp){
            err = tmp;
	}
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
    
    size_t n = A->n;
    // my RHS and the reference one
    double* b1 = calloc(n, sizeof(double));
    double* b2 = calloc(n, sizeof(double));
    
    mesh_RHS(H, b1, F_vol, g_Neu);
    mesh_buildRhs(H, b2, F_vol, g_Neu);

    // check matrix and RHS
    double err_A = vec1_norm(A->x, B->x, A->nzmax);
    printf("Error of matrix is %10g\n", err_A);

    double err_b = vec1_norm(b1, b2, n);
    printf("Error of RHS is %10g\n", err_b);
    
    // TODO: Solve LSE
}
