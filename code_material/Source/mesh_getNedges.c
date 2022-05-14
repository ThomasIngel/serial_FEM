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

  
    /* Set mesh sizes */
  C = In->coord; E = In->elem; B = In->bdry;
  nCoord = In->ncoord; nElem  = In->nelem; nBdry  = In->nbdry;     
  nEdges = 0;                                 /* # of edges */
  for ( k = 0 ; k < nElem ; k++) 
    for (j = 3 ; j < 6 ; j++)
       if ( E[7*k+j] > nEdges) nEdges = E[7*k+j];
  nEdges++;
    
  /* Allocate storage for refined mesh */
  Out = mesh_alloc (nCoord+nEdges, 4*nElem, 2*nBdry);
  In->nedges = nEdges;
  In->edge2no = malloc(2 * nEdges * sizeof(index)); 
  
  edge2no = In->edge2no;
  nC = Out->coord; nE = Out->elem; nB = Out->bdry;
    
  /* Get endpoints for each edge, i.e. compute edgeno */              
  for ( k = 0 ; k < nElem ; k++){
    for ( j = 0 ; j < 3 ; j++){
      p = 7 * k;
      edge2no[2 * E[p+j+3]  ] = E[p + j]; 
      edge2no[2 * E[p+j+3]+1] = E[p + isucc[j]];
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