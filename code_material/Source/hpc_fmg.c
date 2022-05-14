#include "hpc.h"

void hpc_fmg(sed **A, double *b, double *x, index nCycle,
             mesh **H, index nLevel, index pre, index post, index gamma)          
{    
	index j, k, nIter, *p;
  double *wx, *wb, *wr, **xx, **bb, **rr, tmp;  
  
  if( !nLevel ) /* If no hierachy -> solve exactly */
  { 
    wr = malloc( A[0]->n * sizeof(double));
    for (k = 0; k < 5 * A[0]->n; k++){
      sed_gs_constr(A[0], b, x, wr, H[0]->fixed, H[0]->nfixed, 1);
      sed_gs_constr(A[0], b, x, wr, H[0]->fixed, H[0]->nfixed, 0);
    }
    free(wr);
    nIter = 2 * k;
	}
  else
  {
    p = malloc ( (nLevel+2) * sizeof(index)); 
    p[0] = 0;
    for (j = 0; j <=nLevel; j++) p[j+1] = p[j] + A[j]->n;
  
    wx = calloc(p[nLevel], sizeof(double));
    wb = malloc(p[nLevel]   * sizeof(double));
    wr = malloc(p[nLevel+1] * sizeof(double));
  
    xx = malloc( (nLevel+1) * sizeof(double*));
    bb = malloc( (nLevel+1) * sizeof(double*));
    rr = malloc( (nLevel+1) * sizeof(double*));
    
    for (j = 0; j < nLevel; j++){
      xx[j] = wx + p[j]; bb[j] = wb + p[j]; rr[j] = wr + p[j];
    }
    xx[nLevel] = x; bb[nLevel] = b; rr[nLevel] = wr + p[nLevel];
    free(p); 
    
    /* restriction of rhs downwards */
    for ( k = nLevel; k > 0; k--) {
      hpc_rest(bb[k], H[k]->edge2no, H[k]->nedges, bb[k-1], H[k]->ncoord);
    }
    
    // "Solve on coarsest grid"
    for (k = 0; k < 5 * A[0]->n; k++){
      sed_gs_constr(A[0], bb[0], xx[0], rr[0], H[0]->fixed, H[0]->nfixed, 1);
      sed_gs_constr(A[0], bb[0], xx[0], rr[0], H[0]->fixed, H[0]->nfixed, 0);
    }
    
    /* iterate until convergence */
    for ( k = 0 ; k < nLevel; k++) 
    {
      /* prolongation */ 
      hpc_prol_quad(xx[k], xx[k+1], H[k]->elem, H[k]->ncoord, 
                    H[k]->nelem, H[k]->nedges);
      /* do multigrid loops */
      for ( j = 0 ; j < nCycle; j++) 
      hpc_mg_cycle(A, H, k+1, bb, xx, rr, pre, post, gamma);
    }
    free(wx); free(wb); free(wr);
    free(xx); free(bb); free(rr);
	}
	return ;
}

