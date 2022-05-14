#include "hpc.h"
/* Collect all nodes of typ 0 in fixedNodes */
index mesh_getEdge2no(const index nElem, const index *Elem, 
                      index *nEdges, index **edge2no)
{
  index k, j, nE;
  static index isucc[3] = {1,2,0};

  if ( !nEdges ) return(0);
  /* Set mesh sizes */
  nE = 0;                                 /* # of edges */
  for ( k = 0 ; k < nElem ; k++) 
    for (j = 3 ; j < 6 ; j++)
       if ( Elem[7*k+j] > nE) nE = Elem[7*k+j];
  nE++;
  *nEdges = nE;  
  /* Allocate storage for refined mesh */
  *edge2no = malloc(2 * nE * sizeof(index)); 
  if ( !(*edge2no) ) return(0);

  /* Get endpoints for each edge, i.e. compute edgeno */              
  for ( k = 0 ; k < nElem ; k++){
    for ( j = 0 ; j < 3 ; j++){
       (*edge2no)[2 * Elem[7 * k + j + 3]  ] = Elem[7 * k + j]; 
       (*edge2no)[2 * Elem[7 * k + j + 3]+1] = Elem[7 * k + isucc[j]];
    }
  }
  return (1);
}