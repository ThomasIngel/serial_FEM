#include "hpc.h"
/* print a sparse matrix; use %g for integers to avoid differences with index */
index cs_print (const cs *A, index brief)
{
    index p, j, m, n, nzmax, nz, *Ap, *Aind ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Aind = A->ind ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    if (HPC_CSC (A))
    {
        printf ("%g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n])) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    col %g : locations %g to %g\n", (double) j, 
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                printf ("      %g : %g\n", (double) (Aind [p]), Ax ? Ax [p] : 1) ;
                if (brief && p > 10) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else if (HPC_CSR (A))
    {
        printf ("%g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n])) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    row %g : locations %g to %g\n", (double) j, 
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                printf ("      %g : %g\n", (double) (Aind [p]), Ax ? Ax [p] : 1) ;
                if (brief && p > 10) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else 
    {
        printf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            printf ("    %g %g : %g\n", (double) (Aind [p]), (double) (Ap [p]),
                Ax ? Ax [p] : 1) ;
            if (brief && p > 10) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}
