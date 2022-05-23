#include "hpc.h"
#include "mesh_trans.h"
#include <mpi.h>

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

void* make_global(index *c,index nlocal,double *rhs_loc,double *rhs_glob){
  for(int i=0;i<nlocal;i++){
    rhs_glob[c[i]] = rhs_loc[i];
  }
}

void* make_local(index *c,index nlocal,double *rhs_glob,double *rhs_loc){
  for(int i=0;i<nlocal;i++){
    rhs_loc[i] = rhs_glob[c[i]] ;
  }
}

int main(int argc, char *argv[]) {

  int numprocs;
	int myid;
	int i;
  const int root=0;
	MPI_Status stat;

  MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* find out how big the SPMD world is */
	MPI_Comm_rank(MPI_COMM_WORLD,&myid); /* and this processes' rank is */

  mesh_trans **metra;
  index anz_dom = 2;
  index ncoord = 9;
  if (myid == 0){
    
    mesh* H = get_refined_mesh(1);

    double* b1 = calloc(ncoord, sizeof(double));  
    mesh_RHS(H, b1, F_vol, g_Neu); 
    printf("\nProcessor %d rhs full mesh:\n", myid);
    for(int i=0;i<ncoord;i++){
      printf("%lg\n",b1[i]);
    }

    metra = malloc ( (anz_dom) * sizeof(mesh_trans));

    for(size_t i=0;i<anz_dom;i++){
      metra[i]=alloc_mesh_trans(anz_dom,ncoord);
    }

    meshsplit(H, metra, anz_dom);
  }
  mesh_trans* mesh_loc =  scatter_meshes(metra,MPI_COMM_WORLD,anz_dom,ncoord);

  sed* S;
  S = malloc (sizeof(sed));
  S = sed_sm_build(mesh_loc);
  /*printf("\nProcessor %d lokale SM:\n", myid);
  sed_print(S,0);*/
  double* b = calloc(mesh_loc->ncoord_loc, sizeof(double));
  mesh_trans_rhs(mesh_loc,b,F_vol, g_Neu);

  printf("\nProcessor %d rhs:\n", myid);
  for(int i=0;i<mesh_loc->ncoord_loc;i++){
    printf("%lg\n",b[i]);
  }

  double* rhs_global = calloc(mesh_loc->ncoord_glo, sizeof(double));
  make_global(mesh_loc->c,mesh_loc->ncoord_loc,b,rhs_global);

  if(myid==0){
    /*printf("\nProcessor %d rhs:\n", myid);
    for(int i=0;i<mesh_loc->ncoord_loc;i++){
      printf("%lg\n",b[i]);
    }*/
    printf("\nProcessor %d made global rhs:\n", myid);
    for(int i=0;i<mesh_loc->ncoord_glo;i++){
      printf("%lg\n",rhs_global[i]);
    }
  }

  double recv_global[mesh_loc->ncoord_glo];
  MPI_Allreduce(
    rhs_global,
    recv_global,
    9,
    MPI_DOUBLE,
    MPI_SUM,
    MPI_COMM_WORLD);
  double* rhs_local = calloc(mesh_loc->ncoord_loc, sizeof(double));
  make_local(mesh_loc->c,mesh_loc->ncoord_loc,rhs_global,rhs_local);

  if(myid==0){
    printf("\nProcessor %d rhs after allreduce:\n", myid);
    for(int i=0;i<mesh_loc->ncoord_glo;i++){
      printf("%lg\n",recv_global[i]);
    }
    /*printf("\nProcessor %d rhs after local:\n", myid);
    for(int i=0;i<mesh_loc->ncoord_loc;i++){
      printf("%lg\n",rhs_local[i]);
    }*/
  }

  MPI_Finalize();
  return 0;
}
  


