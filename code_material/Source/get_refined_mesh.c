#include "hpc.h"

// Load Mesh from file and refine norefine times
mesh *get_refined_mesh(int norefine){

    // Load and initialize Mesh
    mesh **M;
    M = malloc ( (norefine+1) * sizeof(mesh));
    char *fname = "../code_material/Problem/problem1";
    M[0] = mesh_load(fname);   
    mesh_getEdge2no(M[0]->nelem, M[0]->elem, &M[0]->nedges, &M[0]->edge2no);
    M[0]->fixed = mesh_getFixed(M[0]->ncoord, M[0]->bdry, M[0]->nbdry, &M[0]->nfixed);
    
    // Refine Mesh
    for (int k = 0; k<norefine; k++){
        M[k+1] = mesh_refine(M[k]);
        mesh_getEdge2no(M[k+1]->nelem, M[k+1]->elem, &M[k+1]->nedges, &M[k+1]->edge2no);
        M[k+1]->fixed = mesh_getFixed(M[k+1]->ncoord, M[k+1]->bdry, M[k+1]->nbdry, &M[k+1]->nfixed);
        free(M[k]);
    }
    return(M[norefine]);
}
