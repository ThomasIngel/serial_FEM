#include "hpc.h"
/* print a sparse matrix; */
index mesh_print (const mesh *M, index brief)
{
    index j, k, ncoord, nelem, nbdry, nfixed, *Elem, *Bdry, *Fixed, total ;
    double *Coord ;
    if (!M) { printf ("(null)\n") ; return (0) ; }
    ncoord = M->ncoord ; nelem = M->nelem ; nbdry = M->nbdry ; nfixed = M->nfixed ; 
    Coord = M->coord ; Elem = M->elem ; Bdry = M->bdry ;  Fixed = M->fixed ; 
    printf ("no. of coordinates %zu\nno. of elements: %zu\nno. of bdry elements: %zu\n", 
            ncoord, nelem, nbdry) ;
    printf ("(x,y) - coordinates\n"); 
    for (j = 0 ; j < ncoord ; j++)
    {
      printf ("    (%lg,  %lg)\n", Coord[2*j], Coord[2*j+1] );
      if (brief && j > 10) { printf ("  ...\n") ; break ; }
    }
    printf ("vertices , midpoints, type\n"); 
    for (j = 0 ; j < nelem ; j++)
    {
      for (k = 0 ; k < 7 ; k++) printf (" %zu",Elem[7*j+k]);
      printf ("\n");
      if (brief && j > 10) { printf ("  ...\n") ; break ; }
    }
    printf ("boundray elements (endpoints, midpoints, type)\n"); 
    for (j = 0 ; j < nbdry ; j++)
    {
      for (k = 0 ; k < 4 ; k++) printf (" %zu", Bdry[4*j+k]);
      printf ("\n");
      if (brief && j > 10) { printf ("  ...\n") ; break ; }
    }
    if (nfixed)
    {
      printf ("fixed nodes\n"); 
      for (j = 0 ; j < nfixed ; j++)
      {
        printf (" %zu\n", Fixed[j]);
        if (brief && j > 10) { printf ("  ...\n") ; break ; }
      }
    }
    
    printf ("\nMemory\n");
    printf ("Coordinates : %12zu Byte\n", ncoord*2*sizeof(double));
    printf ("Elements :    %12zu Byte\n", nelem*7*sizeof(index));
    printf ("Boundary :    %12zu Byte\n", nbdry*4*sizeof(index));
    printf ("Edge2no :     %12zu Byte\n", M->nedges*2*sizeof(index));
    total = ncoord*2*sizeof(double) 
          + (7*nelem+4*nbdry+M->nedges*2)*sizeof(index);
    printf ("Total :       %12.6g MByte\n", (double) total/1024./1024.);

    return (1) ;
}
