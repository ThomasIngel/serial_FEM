#include "hpc.h"

// Load Mesh from file and refine norefine times
mesh *get_refined_mesh(int norefine){

    // Load and initialize Mesh
    mesh **M;
    M = malloc ( (norefine+1) * sizeof(mesh));
    char *fname = "../Problem/problem1";
    
    printf("load file from mesh \n");
    M[0] = mesh_load(fname);
    printf("finished loading mesh from file \n");
    printf("----------------------------\n");
    printf("Test\n");
    printf("M = %p \n", (void *) &M[0]);
    /*printf("nelem = %d \n", M[0]->nelem);
    printf("ncoords = %d \n", M[0]->ncoord);      
    mesh_getEdge2no(M[0]->nelem, M[0]->elem, &M[0]->nedges, &M[0]->edge2no);
    M[0]->fixed = mesh_getFixed(M[0]->ncoord, M[0]->bdry, M[0]->nbdry, &M[0]->nfixed);
    */
    // Refine Mesh
    printf("entering loop \n");
    for (int k = 0; k<norefine; k++){
    	printf("starting %d refinement \n", k+1);
        M[k+1] = mesh_refine(M[k]);
        printf("finished mesh_refinement method\n");
        mesh_getEdge2no(M[k+1]->nelem, M[k+1]->elem, &M[k+1]->nedges, &M[k+1]->edge2no);
        M[k+1]->fixed = mesh_getFixed(M[k+1]->ncoord, M[k+1]->bdry, M[k+1]->nbdry, &M[k+1]->nfixed);
        printf("finished the other strange methods\n");
        free(M[k]);
        printf("freed memory from the refinement bevor\n");
    }
    return(M[norefine]);
}
