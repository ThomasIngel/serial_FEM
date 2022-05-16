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

int main() {
	// inital mesh
    mesh* H = get_refined_mesh(1);
    printf("loading done");
    sed* A = sed_nz_pattern(H);
    sed* B = sed_nz_pattern(H);
    
    // construct my stiffnes matrix and print it
    mesh_stima_global(H, A);
    sed_print(A, 0);
    
    // construct stiffnes matrix from the material and print it
    sed_buildS(H, B);
    sed_print(B, 0);
    
    // TODO: Construct RHS
    
    // TODO: Solve LSE
}
