// Adapted sed_spmv - 17.05.22 (mit alpha)

#include "hpc.h"
/* y = alpha * A * x + y */
index sed_spmv_adapt (const sed *A, const double *x, double *y, double alpha)
{
  index p, j, m, n, nz, *Ap, *Ai ;
  double *Ax, tmp ;

  if (!A || !x || !y) return (0) ;                /* check inputs */
  n = A->n ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
  {
    y[j] += alpha * Ax[j] * x[j] ;
    for (p = Ai[j] ; p < Ai[j+1] ; p++)
    {
      y[Ai[p]] += alpha * Ax[p] * x[j];
      y[j] += alpha * Ax[p] * x[Ai[p]];
    }
  }
  return (1) ;
}
