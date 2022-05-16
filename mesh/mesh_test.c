#include "hpc.h"

int main()
{
    mesh **M;
    M = malloc (sizeof(mesh));

    printf("\n========================================\n");

    int norefine = 1;
    M[0] = get_refined_mesh(norefine);
    mesh_print(M[0],1);
    printf("\n========================================\n");

    free(M); 

    return (0) ;
}
