#include "hpc.h"
/* allocate a sparse matrix (triplet form or compressed-col/row form) */
cs *cs_alloc (index m, index n, index nzmax, index values, index typ)
{
    cs *A = calloc (1, sizeof (cs)) ;      /* allocate the cs struct */
    if (!A) return (NULL) ;                /* out of memory */
    A->m = m ;                             /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = HPC_MAX (nzmax, 1) ;
    A->nz = -typ;                          /* allocate triplet or comp.col/row */
    A->p = malloc ( (!typ ? nzmax : n+1) * sizeof (index)) ;
    A->ind = malloc (nzmax * sizeof (index)) ;
    A->x = values ? malloc (nzmax * sizeof (double)) : NULL ;
    return ((!A->p || !A->ind || (values && !A->x)) ? cs_free (A) : A) ;
}

/* change the max # of entries sparse matrix */
index cs_realloc (cs *A, index nzmax)
{
    index ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (HPC_TRIPLET (A)) ? A->nz : (A->p [A->n]);
    nzmax = HPC_MAX (nzmax, 1) ;
    A->ind = hpc_realloc (A->ind, nzmax, sizeof (index), &oki) ;
    if (HPC_TRIPLET (A)) A->p = hpc_realloc (A->p, nzmax, sizeof (index), &okj) ;
    if (A->x) A->x = hpc_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* free a sparse matrix */
cs *cs_free (cs *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->p) ;      /* free the crs struct and return NULL */
    free (A->ind) ;
    free (A->x) ;
    free (A);
    return (NULL) ; 
}

/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, index ok)
{
    free (w) ;                       /* free workspace */
    free (x) ;
    return (ok ? C : cs_free (C)) ;   /* return result if OK, else free it */
}


