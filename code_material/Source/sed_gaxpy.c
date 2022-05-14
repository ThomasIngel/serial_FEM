#include "hpc.h"
/* y = A * x + y */
index sed_gaxpy (const sed *A, const double *x, double *y)
{
  index p, j, m, n, nz, *Ap, *Ai ;
  double *Ax, tmp ;
  
  if (!A || !x || !y) return (0) ;                /* check inputs */
  n = A->n ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
  {
    y[j] += Ax[j] * x[j] ;
    for (p = Ai[j] ; p < Ai[j+1] ; p++)
    {
      y[Ai[p]] += Ax[p] * x[j] ;
      y[j]     += Ax[p] * x[Ai[p]] ;
    }
  }
  return (1) ;
}
