#include "hpc.h"
/* print a sparse matrix; use %g for integers to avoid differences with index */
index sed_print (const sed *A, index brief)
{
    index p, j, n, nzmax, *Ai ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    n = A->n ; Ai = A->i ; Ax = A->x ; nzmax = A->nzmax ;
    
    printf ("%g-by-%g, nzmax: %g nnz: %g\n", (double) n,
            (double) n, (double) nzmax, (double) Ai[n]) ;
    printf ("diagonal entries\n"); 
    for (j = 0 ; j < n ; j++)
    {
       printf ("      %g : %g\n", (double) j, Ax ? Ax [j] : 1) ;
       if (brief && j > 10) { printf ("  ...\n") ; break ; }
    }
      
    for (j = 0 ; j < n ; j++)
    {
      printf ("    col %g : locations %g to %g\n", (double) j, 
                (double) (Ai [j]), (double) (Ai [j+1]-1)) ;
      for (p = Ai [j] ; p < Ai [j+1] ; p++)
      {
        printf ("      %g : %g\n", (double) (Ai [p]), Ax ? Ax [p] : 1) ;
        if (brief && p > 10) { printf ("  ...\n") ; return (1) ; }
      }
    }
    return (1) ;
}
