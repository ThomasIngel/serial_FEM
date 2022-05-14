#include "hpc.h"
/* Collect all nodes of typ 0 in fixedNodes */
index *mesh_getFixed(const index nCoord, const index *bdry, 
                     const index nBdry, index *nFixed)
{
  index i, j, cnt, nz, *fixed, ok = 1, MAX_IND = 0;
  bool *flag;

  for ( i=0 ; i < nBdry ; i++)
  {
    if (!bdry[4*i+3]){
      MAX_IND = HPC_MAX(MAX_IND, bdry[4*i+2]);
    }
  }
  MAX_IND += nCoord + 1;
  flag = (bool*) calloc(MAX_IND,sizeof(bool));
  cnt = 0;
  for ( i=0 ; i < nBdry ; i++){
    if (!bdry[4*i+3]){
      cnt++;
      flag[bdry[4*i]]   = 1;
      flag[bdry[4*i+1]] = 1;
      flag[nCoord+bdry[4*i+2]] = 1;
    }
  }
  fixed = malloc(3*cnt * sizeof(index));
  if (!fixed) return (NULL);
  nz = 0;
  for ( j = 0 ; j < MAX_IND; j++) 
  {
    if (flag[j]) fixed[nz++] = j;
  }
  fixed = realloc(fixed,nz*sizeof(index));
  *nFixed = nz;
  free(flag);
  
  return (fixed);
}