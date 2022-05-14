#include "hpc.h"

index hpc_mg(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H, index nLevel, index pre, index post, index gamma)          
{
    
  index j, k, nIter, *p;
  double *wx, *wb, *wr, **xx, **bb, **rr, tmp;
  
  if( !nLevel ) /* If no hierachy -> solve exactly */
  { 
//	  nIter = hpc_cg_constr(A[0], b, x, H[0]->fixed, 
//                          H[0]->nfixed, tol, maxit);
    wr = malloc( A[0]->n * sizeof(double));
    for (k = 0; k < 10; k++){
      sed_gs_constr(A[0], b, x, wr, H[0]->fixed, H[0]->nfixed, 1);
      sed_gs_constr(A[0], b, x, wr, H[0]->fixed, H[0]->nfixed, 0);
    }
    free(wr);
    nIter = 2 * k;
	}
  else
  {
    p = malloc ((nLevel+2) * sizeof(index)); 
    p[0] = 0;
    for (j = 0; j <=nLevel; j++) p[j+1] = p[j] + A[j]->n;
  
    wx = malloc( p[nLevel]   * sizeof(double) );
    wb = malloc( p[nLevel]   * sizeof(double) );
    wr = malloc( p[nLevel+1] * sizeof(double) );
  
    xx = malloc( (nLevel+1) * sizeof(double*) );
    bb = malloc( (nLevel+1) * sizeof(double*) );
    rr = malloc( (nLevel+1) * sizeof(double*) );
    
    for (j = 0; j < nLevel; j++){
      xx[j] = wx + p[j]; bb[j] = wb + p[j]; rr[j] = wr + p[j];
    }
    xx[nLevel] = x; bb[nLevel] = b; rr[nLevel] = wr + p[nLevel];
    free(p); 
    
    /* iterate until convergence */
    for ( j = 0 ; j <= maxit; j++) 
    {
      /* Compute residual, r = b - A * x */
      for ( k = 0 ; k < A[nLevel]->n ; k++) rr[nLevel][k] = 0.0;
      sed_gaxpy(A[nLevel], xx[nLevel], rr[nLevel]);
      for ( k = 0 ; k < A[nLevel]->n ; k++) rr[nLevel][k] = bb[nLevel][k] - rr[nLevel][k];
      /* Consider constrains */
      for ( k = 0 ; k < H[nLevel]->nfixed; k++) 
            rr[nLevel][H[nLevel]->fixed[k]] = 0.0;  
        
      tmp = 0.0;
      for ( k=0 ; k < A[nLevel]->n ; k++) tmp += rr[nLevel][k] * rr[nLevel][k];
               
      if ( tmp <= tol*tol ) 
      {
        free(wx); free(wb); free(wr);
        free(xx); free(bb); free(rr);
        return (j);
      }
      hpc_mg_cycle(A, H, nLevel, bb, xx, rr, pre, post, gamma);
    }
    fprintf(stderr,
        "Max iterations reached: maxit = %ld, r = %g\n", maxit, sqrt(tmp));
    nIter = j;
    free(wx); free(wb); free(wr);
    free(xx); free(bb); free(rr);
	}
	return (nIter);
}