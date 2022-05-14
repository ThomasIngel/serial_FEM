#include "hpc.h"
/* cg-olver with constrains */
void hpc_prol(double *x, index nx, index *edgeno, index nEdges, double *y)
{ 
  index i;

  for (i = 0; i < nx; i++) y[i] += x[i];  
  for (i = 0; i < nEdges; i++){
    y[nx+i] += 0.5 * (x[edgeno[2*i  ]]
                    + x[edgeno[2*i+1]]);
  }
  return;
}
                                 