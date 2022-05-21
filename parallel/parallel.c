#include "hpc.h"
#include "mesh_trans.h"
#include <mpi.h>

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
  printf("\nProcessor %d lokale SM:\n", myid);
  sed_print(S,0);

  MPI_Finalize();
  return 0;
}
  


