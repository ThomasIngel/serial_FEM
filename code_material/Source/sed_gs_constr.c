#include "hpc.h"

/* Gauss-Seidel Iteration with constrains x(fixed) = b(fixed) */
index sed_gs_constr (const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed, index forward)
{
    index p, j, m, n, nz, t, ft, *Ap, *Ai ;
    double *Ax, sum, den ;
    
    if ( !A || !x || !b ) return (0) ;  /* check inputs */
    n = A->n ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n; j++) w[j] = b[j];
    
    if (forward) 
    {
      for (j = 0 ; j < n; j++)
      {
        for (p = Ai[j] ; p < Ai[j+1] ; p++) w[j] -= Ax[p] * x [Ai[p]];
      }
      ft = fixed[0]; t = 1; 
      for (j = 0 ; j < n; j++)
      {
        if ( j != ft)
        { 
          x[j] = w[j] / Ax[j];
        } 
        else 
        {
          if (t < nFixed) ft = fixed[t++];
        }
        for (p = Ai[j] ; p < Ai[j+1] ; p++) w[Ai[p]] -= Ax[p] * x [j];    
      }
    }
    else
    {
      for (j = 0 ; j < n; j++)
      {
        for (p = Ai[j] ; p < Ai[j+1] ; p++) w[Ai[p]] -= Ax[p] * x [j];
      }
      t = nFixed-1; ft = fixed[nFixed-1]; 
      for (j = n-1 ; j >=0; j--)
      {
        if ( j != ft)
        {  
          for (p = Ai[j] ; p < Ai[j+1] ; p++)
          {
           w[j] -= Ax[p] * x [Ai[p]];
          }
          x[j] = w[j] / Ax[j];
        } 
        else 
        {
          if (t > 0 ) ft = fixed[--t];
        } 
      }
    }    
//    for (j = 0 ; j < nFixed; j++) x[fixed[j]] = b[fixed[j]];
    return (1) ;
}
