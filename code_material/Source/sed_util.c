#include "hpc.h"
/* allocate a sparse square matrix with extracted diagonal in compressed col. form) */
sed *sed_alloc (index n, index nzmax, index values)
{
    sed *A = calloc (1, sizeof (sed)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                /* out of memory */
    A->n = n ;                             /* define dimensions and nzmax */
    A->nzmax = nzmax = HPC_MAX (nzmax, n+1) ;
    A->i = malloc (nzmax * sizeof (index)) ;
    if (!A->i) return (NULL);
    A->i[A->n] = nzmax;
    A->x = values ? malloc (nzmax * sizeof (double)) : NULL ;
    return ((!A->i || (values && !A->x)) ? sed_free (A) : A) ;
}

/* change the max # of entries sparse matrix */
index sed_realloc (sed *A, index nzmax)
{
    index ok, oki, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = A->i[A->n];
    nzmax = HPC_MAX (nzmax, 1+A->n) ;
    A->i = hpc_realloc (A->i, nzmax, sizeof (index), &oki) ;
    if (A->x) A->x = hpc_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* free a sparse matrix */
sed *sed_free (sed *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->i) ;      /* free the sed struct and return NULL */
    free (A->x) ;
    free (A);
    return (NULL) ; 
}

/* free workspace and return a sparse matrix result */
sed *sed_done (sed *C, void *w, void *x, index ok)
{
    free (w) ;                       /* free workspace */
    free (x) ;
    return (ok ? C : sed_free (C)) ;   /* return result if OK, else free it */
}


