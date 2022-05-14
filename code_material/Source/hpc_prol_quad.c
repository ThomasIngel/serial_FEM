#include "hpc.h"
/* cg-olver with constrains */
void hpc_prol_quad(double *x, double *y, index *elem, index nC, index nT, index nE)
{ 
  index ia, ib, ic, ma, mb, mc, j, k, isucc[3] = {1,2,0}, iprae[3] = {2,0,1};
  
  for (j = 0; j < nC + nE; j++)
  {
    y[j] = x[j] ;
  }
  for (j = 0; j < nT; j++){
    for (k = 0; k < 3; k++){
      ia = elem[7*j+k];   ib = elem[7*j+isucc[k]];   ic = elem[7*j+iprae[k]];
      ma = elem[7*j+k+3]; mb = elem[7*j+isucc[k]+3]; mc = elem[7*j+iprae[k]+3];
      if ( ia < ib )
      {
        y[nC+nE+2*ma]   = (3.0 * x[ia] - x[ib] + 6.0 * x[nC+ma]) / 8.0;
        y[nC+nE+2*ma+1] = (3.0 * x[ib] - x[ia] + 6.0 * x[nC+ma]) / 8.0;
      }
      else
      {
        y[nC+nE+2*ma+1] = (3.0 * x[ia] - x[ib] + 6.0 * x[nC+ma]) / 8.0;
        y[nC+nE+2*ma  ] = (3.0 * x[ib] - x[ia] + 6.0 * x[nC+ma]) / 8.0;
      }
      y[nC+3*nE+3*j+k] = (4*x[nC+ma]+2*x[nC+mb]+4*x[nC+mc]-x[ib]-x[ic])/8.0; 
    }
  }
  return;
}
                                 