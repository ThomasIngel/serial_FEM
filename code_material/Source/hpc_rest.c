#include "hpc.h"
/* cg-olver with constrains */

void hpc_rest(double *x, index *edgeno, index nEdges, double *y, index ny)
{ 
  index i;
  
  for (i=0; i<ny; i++) y[i] = x[i];  
  for (i=0; i<nEdges; i++){
    y[edgeno[2*i  ]] += 0.5*x[ny+i];
    y[edgeno[2*i+1]] += 0.5*x[ny+i];
  }
  return;
}

