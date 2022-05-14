#include "hpc.h"

/* wrapper for free */
void *hpc_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of hpc_free */
}

/* wrapper for realloc */
void *hpc_realloc (void *p, index n, size_t size, index *ok)
{
    void *pnew ;
    pnew = realloc (p, HPC_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}
