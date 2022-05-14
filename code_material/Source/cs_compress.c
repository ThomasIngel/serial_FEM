#include "hpc.h"
/* C = compressed-column form of a triplet matrix T */
cs *cs_compress (const cs *T, index typ) /* typ 1 -> csc, typ 2 -> csr */
{
    index m, n, sizep, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    double *Cx, *Tx ;
    cs *C ;
    if (!HPC_TRIPLET (T)) return (NULL) ;                /* check inputs */
    
    m = T->m ; n = T->n; Tx = T->x ; nz = T->nz ;
    if (typ == 1)
    {
      Ti = T->ind; Tj = T->p; sizep = n ;
    }
    else 
    {
      Ti = T->p; Tj = T->ind; sizep = m ;
    }
    C = cs_alloc (m, n, nz, Tx != NULL, typ) ;       /* allocate result */
    w = calloc (sizep, sizeof (index)) ;            /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;   /* out of memory */    
    Cp = C->p ; Ci = C->ind ; Cx = C->x ;
    for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;           /* column counts */
    hpc_cumsum (Cp, w, sizep) ;                         /* column pointers */
    for (k = 0 ; k < nz ; k++)
    {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }
    return (cs_done (C, w, NULL, 1)) ;     /* success; free w and return C */
}
