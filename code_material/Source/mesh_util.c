#include "hpc.h"

/* allocate mesh data structure */
mesh *mesh_alloc (index ncoord, index nelem, index nbdry)
{
    mesh *M = malloc (sizeof (mesh)) ;  /* allocate the mesh struct */
    if (!M) return (NULL) ;                    /* out of memory */
    M->ncoord = ncoord ;                        
    M->nelem = nelem ;
    M->nbdry = nbdry;
    M->nedges = 0;
    M->nfixed = 0;
    M->coord = malloc (ncoord * 2 * sizeof (double)) ;
    M->elem = malloc (nelem * 7 * sizeof (index)) ;
    M->bdry = malloc (nbdry * 4 * sizeof (index)) ;
    M->edge2no = NULL; M->fixed=NULL;
    return ((!M->coord || !M->elem || !M->bdry) ? mesh_free (M) : M) ;
}

/* free a mesh data structure */
mesh *mesh_free (mesh *M)
{
    if (!M) return (NULL) ;     /* do nothing if M already NULL */
    free (M->coord) ;
    free (M->elem) ;
    free (M->bdry) ;
    if (M->edge2no) free (M->edge2no) ;
    if (M->fixed)   free (M->fixed) ;
    free (M);
    return (NULL) ;   /* free the mesh struct and return NULL */
}