#include "hpc.h"

// Load Hesh from file and refine norefine times
mesh *get_refined_mesh(int norefine){

    // Load and initialize Hesh
    mesh* I;
    mesh* H;

    I = malloc (sizeof(mesh));
    char *fname = "../code_material/Iroblem/problem1";
    I = mesh_load(fname);   
    mesh_getEdge2no(I->nelem, I->elem, &I->nedges, &I->edge2no);
    I->fixed = mesh_getFixed(I->ncoord, I->bdry, I->nbdry, &I->nfixed);
    
    // Refine Hesh
    for (int k = 0; k<norefine; k++){
        H = mesh_refine(I);
        mesh_getEdge2no(H->nelem, H->elem, &H->nedges, &H->edge2no);
        H->fixed = mesh_getFixed(H->ncoord, H->bdry, H->nbdry, &H->nfixed);
        I = mesh_free(I);
        I = H;
    }
    return(H);
}
