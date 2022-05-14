#include "hpc.h"
/* remove duplicate entries from A */
index sed_dupl (sed *A)
{
    index i, j, p, q, nz = 0, n, m, *Ai, *w ;
    double *Ax ;
    if (!A) return (0) ;           /* check inputs */
    n = A->n ; Ai = A->i ; Ax = A->x ;
    w = malloc (n * sizeof (index)) ;           /* get workspace */
    if (!w) return (0) ;                        /* out of memory */
    for (i = 0 ; i < n ; i++) w [i] = n ;       /* row i not yet seen */
    nz = Ai[0];
    for (j = 0 ; j < n ; j++)
    {
        q = nz ;                                /* column j will start at q */
        for (p = Ai [j] ; p < Ai [j+1] ; p++)
        {
            i = Ai [p] ;                        /* A(i,j) is nonzero */
            if (w [i] >= q)
            {
                if (Ax) Ax[w[i]] += Ax [p] ;    /* A(i,j) is a duplicate */
            }
            else
            {
                w [i] = nz ;                    /* record where row i occurs */
                Ai [nz] = i ;                   /* keep A(i,j) */
                if (Ax) Ax[nz] = Ax [p] ;
                nz++;
            }
        }
        Ai [j] = q ;                            /* record start of column j */
    }
    Ai [n] = nz ;                               /* finalize A */
    free (w) ;                                  /* free workspace */
    return (sed_realloc(A,-1)) ;                   /* remove extra space from A */
}
