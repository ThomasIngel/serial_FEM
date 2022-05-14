#ifndef _HPC_H
#define _HPC_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>

#include <stdbool.h>

#define index ptrdiff_t

/* --- primary HPC routines and data structures ------------------------- */


typedef struct sed_sparse  /* matrix in sparse matrix in compressed col. */
{                          /* with extracted diagonal storage form      */
    index nzmax ;     /* maximum number of entries */
    index   n ;       /* number of rows/columns          */
    index  *i ;       /* col pointers and row indices    */
    double *x ;       /* numerical values, size i[n] */
} sed ;

typedef struct mesh_data  /* mesh */
{
    index ncoord ;    /* number of coordinates  */
    index nelem ;     /* number of elements   */
    index nedges ;    /* number of edges  */
    index nbdry ;     /* number of boundary elements  */
    index nfixed;     /* number of fixed nodes ????    */
    double *coord ;   /* coordinates (x1,y1,x2,y2, ... ,x_ncoord,y_ncoord) */
    index *elem ;     /* elements ([e1,e2,e3,m1,m2,m3,t1], ... ) */
    index *edge2no ;  /*  */
    index *bdry ;     /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
    index *fixed ;    /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
} mesh ;

/* utilities */
void *hpc_realloc (void *p, index n, size_t size, index *ok);
double hpc_cumsum (index *p, index *c, index n);
 
sed *sed_alloc (index n, index nzmax, index values);
index sed_realloc (sed *A, index nzmax);
sed *sed_free (sed *A);
sed *sed_done (sed *C, void *w, void *x, index ok);
// sed *sed_compress (const cs *A);
index sed_print (const sed *A, index brief);
index sed_gaxpy (const sed *A, const double *x, double *y);
index sed_dupl (sed *A);
index sed_gs_constr (const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed, index forward);

mesh *mesh_alloc (index ncoord, index nelem, index nbdry);
mesh *mesh_free (mesh *M);
mesh *mesh_load (char *fname);
index *mesh_getFixed(const index nCoord, const index *bdry, 
                     const index nBdry, index *nFixed);
index mesh_print (const mesh *M, index brief);
mesh *mesh_refine(mesh *In);
index mesh_getEdge2no(const index nElem, const index *Elem, 
                      index *nEdges, index **edge2no);




void stima_laplace3(double p1[2], double p2[2], double p3[2],
                    index  typ, double dx[6], double ax[9]);

sed *sed_nz_pattern(mesh *M) ; 
index sed_buildS(mesh *M, sed *T);
void mesh_buildRhs(const mesh *M, double *b, double (*f)(double *, index), 
                   double (*g)(double *, index));

void hpc_fmg(sed **A, double *b, double *x, index nCycle,
             mesh **H, index nLevel, index pre, index post, index gamma);          
index hpc_mg(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H, index nLevel, index pre, index post, index gamma);          
index hpc_mg_cycle(sed **A, mesh **H, index nLevel, 
                   double **b, double **x, double **r,
                   index pre, index post, index gamma);

void hpc_rest(double *x, index *edgeno, index nEdges, double *y, index ny);
void hpc_prol(double *x, index nx, index *edgeno, index nEdges, double *y);
void hpc_prol_quad(double *x, double *y, index *elem, index nC, index nT, index nE);


double kappa( double x[2], index typ );
double F_vol( double x[2], index typ );

#define HPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))
#endif
