#include "hpc.h"

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
	mesh *H;
	mesh *H1;
	sed *A;
	char *fname = "../Problem/problem1";
	
	H = mesh_load(fname);
	printf("-----------------------------------------\n");
	mesh_getEdge2no(H->nelem, H->elem, &H->nedges, &H->edge2no);
	H->fixed = mesh_getFixed(H->ncoord, H->bdry, H->nbdry, &H->nfixed);
	H1 = mesh_refine(H);
	mesh_getEdge2no(H1->nelem, H1->elem, &H1->nedges, &H1->edge2no);
	H1->fixed = mesh_getFixed(H1->ncoord, H1->bdry, H1->nbdry, &H1->nfixed);
	
	A = sed_nz_pattern(H1);
	sed_buildS(H1, A);
	
	
	int n = A->n;
    double* x = calloc (n, sizeof(double));       /* get workspace for sol*/
    double* w = calloc (n, sizeof(double));       /* get temporary workspace */
    double* b = calloc (n, sizeof(double));       /* get workspace for rhs*/
    mesh_buildRhs(H1, b, F_vol, g_Neu); /* build rhs (volume and Neumann data */
    
    /* incorporate Dirichlet data */
    index ncoord = H1->ncoord; 
    index nelem = H1->nelem; 
    index nbdry = H1->nbdry; 
    index *bdry = H1->bdry; 
    H1->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H1->nfixed);
    index nfixed = H1->nfixed ; 
    index nedges = H1->nedges; 
    double* Coord = H1->coord;
	
	
	double x1[2];
	double x2[2];
	double x3[2];
	double m[2];
	
    for (int k = 0; k < nbdry; k++)
    {
      if (!bdry[4*k+3])
      {
        x1[0] = Coord[2 * bdry[4*k]];   x1[1] = Coord[2 * bdry[4*k]+1];
        x2[0] = Coord[2 * bdry[4*k+1]]; x2[1] = Coord[2 * bdry[4*k+1]+1];
        m[0] = ( x1[0] + x2[0] ) / 2.0 ;
        m[1] = ( x1[1] + x2[1] ) / 2.0 ;
        x[bdry[4*k]  ] = u_D(x1);
        x[bdry[4*k+1]] = u_D(x2);
        x[ncoord + bdry[4*k+2]] = u_D(m);
      }
    }
	return 0;
}
