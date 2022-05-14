#include "hpc.h"
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
index cs_entry (cs *T, index i, index j, double x)
{
    if (!HPC_TRIPLET (T) || i < 0 || j < 0) return (0) ;     /* check inputs */
    if (T->nz >= T->nzmax && !cs_realloc (T,2*(T->nzmax))) return (0) ;
    if (T->x) T->x [T->nz] = x ;
    T->ind [T->nz] = i ;
    T->p [T->nz++] = j ;
    T->m = HPC_MAX (T->m, i+1) ;
    T->n = HPC_MAX (T->n, j+1) ;
    return (1) ;
}
