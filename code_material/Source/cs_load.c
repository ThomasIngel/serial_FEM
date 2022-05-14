#include "hpc.h"
/* load a triplet matrix from a file */
cs *cs_load (FILE *f, index issym)
{
  double i, j, x;
  index k, offset, idim, jdim;
  cs *T ;
  if (!f) return (NULL) ;                             /* check inputs */
  T = cs_alloc (0, 0, 1, 1, 0) ;                   /* allocate result */
  offset = 1024; idim = 0; jdim = 0;
  
  rewind(f);                                          
  while (fscanf (f, "%lg %lg %lg\n", &i, &j, &x) == 3)
  {
    offset = HPC_MIN(offset,i); idim = HPC_MAX(idim,i);
    offset = HPC_MIN(offset,j); jdim = HPC_MAX(jdim,j);
    if (!cs_entry (T, (index) i, (index) j, x)) return (cs_free (T)) ;
    if ( issym && (i != j) )
    {
      if (!cs_entry (T, (index) j, (index) i, x)) return (cs_free (T)) ;
    }
  }
  if (offset) {
    for (k = 0 ; k < T->nz; k++)
    {
      T->ind [k] -= offset ;
      T->p [k]   -= offset ;
    }
  }
  T->m = idim-offset+1;
  T->n = jdim-offset+1;
  return (T) ;
}
