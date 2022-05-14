#include "hpc.h"

index hpc_mg_cycle(sed **A, mesh **H, index nLevel, 
                   double **b, double **x, double **r,
                   index pre, index post, index gamma)
{  
	index k, nIter;
  double tmp;
	
	if( !nLevel ) /*If coarsest grid is achieved -> solve exactly */
  { 
    for (k = 0; k < 40*A[0]->n; k++){
      sed_gs_constr(A[0], b[0], x[0], r[0], H[0]->fixed, H[0]->nfixed, 1);
      sed_gs_constr(A[0], b[0], x[0], r[0], H[0]->fixed, H[0]->nfixed, 0);
    }
    nIter = 2*k;
	}
  else
  {
    /* Perform some pre - smoothing steps */
    for ( k = 0 ; k < pre ; k++)
       sed_gs_constr(A[nLevel], b[nLevel], x[nLevel], r[nLevel], 
                     H[nLevel]->fixed, H[nLevel]->nfixed,1);
		
    /* Compute residual, r = b - A * x */
    for ( k = 0 ; k < A[nLevel]->n ; k++) r[nLevel][k] = 0.0;
    sed_gaxpy(A[nLevel], x[nLevel], r[nLevel]);
    for ( k = 0 ; k < A[nLevel]->n ; k++) r[nLevel][k] = b[nLevel][k] - r[nLevel][k];
    
    /* Consider constrains */
    for ( k = 0 ; k < H[nLevel]->nfixed; k++) r[nLevel][H[nLevel]->fixed[k]] = 0.0;                 
		
		/* Restrict the residual -> coarser grid on Level+1 */
    hpc_rest(r[nLevel], H[nLevel]->edge2no, H[nLevel]->nedges, b[nLevel-1], 
             H[nLevel]->ncoord);

    /* Init x on coarser level */
    for ( k = 0 ; k < A[nLevel-1]->n ; k++) x[nLevel-1][k] = 0.0;
		
		/* Compute the correction on the coarser grid (Level+1) */
 		for( k = 0; k < (index) gamma; k++){ 
		    nIter = hpc_mg_cycle(A, H, nLevel-1, b, x, r, pre, post, gamma); 
		}
    
		/* Prolongate the correction to the finer grid (Level) */
    hpc_prol(x[nLevel-1], H[nLevel]->ncoord, H[nLevel]->edge2no, 
             H[nLevel]->nedges, x[nLevel]);
     		
    /* Perform some post - smoothing steps */
    for ( k = 0 ; k < pre ; k++)
      sed_gs_constr(A[nLevel], b[nLevel], x[nLevel], r[nLevel], 
                  H[nLevel]->fixed, H[nLevel]->nfixed,0);
	}
	return (nIter);
}

