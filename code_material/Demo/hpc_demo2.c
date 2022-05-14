#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv[100];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

double kappa( double x[2], index typ )
{
  return ( 1.0 );
}

double F_vol( double x[2], index typ )
{
  return ( 0.0 );
}

double g_Neu( double x[2], index typ )
{
  return ( x[0] * x[1] );
}

double u_D( double x[2])
{
//  return ( 0.0 );
  return ( x[0] * x[1] );
}


int main (int argc, char **argv)
{
    index n, k, ncoord, nelem, nbdry, nfixed, nedges, total, *bdry,
          cnt = 0, N = 0, MAX_ITER = 100;
    char *Pdir = "../Problem/", fname[64];
    double *b, *x, *w, *Coord, x1[2], x2[2], m[2];
    mesh **H, *T ;
    sed **A;
    
    TIME_SAVE(0);
    printf("\n========================================\n");
    if (argc < 2 ){ printf("Problem not specified\n"); return(1); } 
    sprintf(fname,"%s%s",Pdir,argv[1]); /* get problem as parameter */
    if (argc>2){                        /* get no. of refinements  */ 
      if ( (atoi(argv[2]) >0) & (atoi(argv[2]) < 13)) N = atoi(argv[2]);
    } 
    printf("Load data form %s, no. refinements = %g\n", fname, (double) N);   
        
    /* Allocate memory for hierachy */
    H = malloc ( (N+1) * sizeof(mesh));
    A = malloc ( (N+1) * sizeof(sed*));
    
    /* Load problem */
    H[0] = mesh_load (fname);              /* load geometry */
    mesh_getEdge2no(H[0]->nelem, H[0]->elem, &H[0]->nedges, &H[0]->edge2no);
    H[0]->fixed = mesh_getFixed(H[0]->ncoord, H[0]->bdry, H[0]->nbdry, &H[0]->nfixed);
    printf("\nInit mesh  # dofs =  %10g\n",(double)  H[0]->ncoord+H[0]->nedges);
    
    /* Build stiffness matrix, refine mesh and create hierachy  */ 
    k = 0;
    while(1)
    {  
    TIME_SAVE(5+5*k);
      A[k] = sed_nz_pattern(H[k]) ;            /* get pattern of matrix */
      if (!A[k]) return(1);
    TIME_SAVE(6+5*k);
      if ( !sed_buildS(H[k],A[k]) ) return(1); /* assemble coefficient matrix */    
    TIME_SAVE(7+5*k);
      if (k >= N) break;
      H[k+1] = mesh_refine(H[k]);
    TIME_SAVE(8+5*k);
  
      mesh_getEdge2no(H[k+1]->nelem, H[k+1]->elem,
                      &H[k+1]->nedges, &H[k+1]->edge2no);
    TIME_SAVE(9+5*k);
      H[k+1]->fixed = mesh_getFixed(H[k+1]->ncoord, H[k+1]->bdry, 
                                    H[k+1]->nbdry, &H[k+1]->nfixed);
      k++;
    }
    TIME_SAVE(1);
    printf("Final mesh # dofs =  %10g\n",(double)  H[N]->ncoord+H[N]->nedges);
    printf("# refinements     =  %10g\n",(double)  N);

    n = A[N]->n;
    x = calloc (n, sizeof(double));       /* get workspace for sol*/
    w = calloc (n, sizeof(double));       /* get temporary workspace */
    b = calloc (n, sizeof(double));       /* get workspace for rhs*/
    mesh_buildRhs(H[N], b, F_vol, g_Neu); /* build rhs (volume and Neumann data */
    TIME_SAVE(2);
    
    /* incorporate Dirichlet data */
    ncoord = H[N]->ncoord ; nelem = H[N]->nelem ; nbdry = H[N]->nbdry ; 
    bdry = H[N]->bdry; H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
    nfixed = H[N]->nfixed ; nedges = H[N]->nedges; Coord = H[N]->coord;

    for ( k = 0; k < nbdry; k++)
    {
      if (!bdry[4*k+3])
      {
        x1[0] = Coord[2 * bdry[4*k]];   x1[1] = Coord[2 * bdry[4*k]+1];
        x2[0] = Coord[2 * bdry[4*k+1]]; x2[1] = Coord[2 * bdry[4*k+1]+1];
        m[0] = ( x1[0] + x2[0] ) / 2.0 ;
        m[1] = ( x1[1] + x2[1] ) / 2.0 ;
        x[bdry[4*k]  ] = u_D(x1);
        x[bdry[4*k+1]] = u_D(x2);
        x[ncoord + bdry[4*k+2]] = u_D(m);
      }
    }
    
    TIME_SAVE(3);
//    cnt = hpc_mg(A, b, x, 1e-10, 50, H, N, 2, 2, 1); 
    hpc_fmg(A, b, x, 1, H, N, 2, 2, 1);          

    TIME_SAVE(4);
                                   
    for (k=0; k<HPC_MIN(10,A[N]->n); k++){ printf(" x[%g] = %g\n",(double) k, x[k]);}
     
    printf("\n");
    printf("Time load & create hierarchy = %9i ns\n", (int) TIME_ELAPSED(0,1));
    printf("Time building rhs            = %9i ns\n", (int) TIME_ELAPSED(1,2));
    printf("Time Dirichlet values        = %9i ns\n", (int) TIME_ELAPSED(2,3));
    printf("Time solve LSE               = %9i ns\n", (int) TIME_ELAPSED(3,4));
    printf("No. iterations               = %9g\n", (double) cnt);
    printf("========================================\n\n");

    k=0;
    while(1){  
      printf("Pattern    = %9i ns\n", (int) TIME_ELAPSED(5+5*k, 6+5*k));
      printf("Build S    = %9i ns\n", (int) TIME_ELAPSED(6+5*k, 7+5*k));
      if (k >= N) break;
      printf("-------------------------------------\n");
      printf("Refine     = %9i ns\n", (int) TIME_ELAPSED(7+5*k, 8+5*k));
      printf("getEdge2no = %9i ns\n", (int) TIME_ELAPSED(8+5*k, 9+5*k));
      printf("getFixed   = %9i ns\n", (int) TIME_ELAPSED(9+5*k,10+5*k));
      k++;
    }
     
    printf ("\nMemory\n");
    printf ("Coordinates : %12zu Byte\n", ncoord*2*sizeof(double));
    printf ("Elements :    %12zu Byte\n", nelem*7*sizeof(index));
    printf ("Boundary :    %12zu Byte\n", nbdry*4*sizeof(index));
    printf ("Edge2no :     %12zu Byte\n", nedges*2*sizeof(index));
    total = ncoord*2*sizeof(double) 
          + (7*nelem+4*nbdry+nedges*2)*sizeof(index);
    printf ("Total :       %12.6g MByte\n", (double) total/1024./1024.);

    for (k=0; k<=N; k++) {mesh_free(H[k]); sed_free(A[k]);}
    free(H); free(A);

    return (0) ;
}
