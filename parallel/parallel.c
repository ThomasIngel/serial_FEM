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
  mesh_trans* test =  scatter_meshes(metra,MPI_COMM_WORLD,anz_dom,ncoord);

  sed* S;
  S = malloc (sizeof(sed));
  S = sed_sm_build(test);
  /*printf("\nProcessor %d lokale SM:\n", myid);
  sed_print(S,0);*/
  double* b = calloc(test->ncoord_loc, sizeof(double));
  mesh_trans_rhs(test,b,F_vol, g_Neu);

  printf("\nProcessor %d rhs:\n", myid);
  for(int i=0;i<test->ncoord_loc;i++){
    printf("%lg\n",b[i]);
  }

  MPI_Finalize();
  return 0;
}
  


