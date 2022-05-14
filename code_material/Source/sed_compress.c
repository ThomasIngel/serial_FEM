#include "hpc.h"
/* C = sparse compressed-column extracted diagonal form 
 *     of a triplet matrix T */
sed *sed_compress (const cs *T) 
{
    index m, n, sizep, nz, nzmax, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    double *Cx, *Tx ;
    sed *C ;
    if (!HPC_TRIPLET (T)) return (NULL) ;              /* check inputs */
    
    m = T->m ; n = T->n ;  Ti = T->ind; Tj = T->p; 
    Tx = T->x ; nz = T->nz ;
    if (m != n) return (NULL) ;
    
    w = calloc (n, sizeof (index)) ;                   /* get workspace */
    if (!w) return (NULL) ;                            /* out of memory */    
    for (k = 0 ; k < nz ; k++)
    {
      if ( Ti [k] != Tj [k] ) w [Tj [k]]++ ;           /* column counts */
    }
    nzmax = n+1;
    for (k = 0 ; k < n ; k++) nzmax+=w[k];

    C = sed_alloc (n, nzmax, Tx != NULL) ;             /* allocate result */
    if (!C) return (sed_done (C, w, NULL, 0)) ;        /* out of memory */    
    Ci = C->i ; Cx = C->x ;

    hpc_cumsum (Ci, w, n) ;                            /* column pointers */
    for (k = 0 ; k < n ; k++) {Ci[k]+=n+1; w[k]+=n+1;}
    Ci[n]+=n+1;
    for (k = 0 ; k < nz ; k++)
    {
      if ( Ti [k] == Tj [k] )
      {
        if (Cx) Cx [Ti [k]] = Tx [k] ;        
      } 
      else 
      {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
      }        
    }
    return (sed_done (C, w, NULL, 1)) ;     /* success; free w and return C */
}
