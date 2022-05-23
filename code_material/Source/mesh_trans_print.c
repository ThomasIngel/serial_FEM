#include "hpc.h"
#include "mesh_trans.h"


void mesh_trans_print (mesh_trans** metra, index domain)
{
printf("\n------ Mesh of domain %td -----\n", domain);
printf("\nncoord_glob (# Koordinaten global) = %td", metra[domain]->ncoord_glo);
printf("\nncoord_loc (# Koordinaten lokal) = %td", metra[domain]->ncoord_loc);
printf("\nnelem_loc (# Elemente lokal) = %td", metra[domain]->nelem_loc);
printf("\nnbdry_loc (# Gebietsrandpunkte im Gebiet) = %td", metra[domain]->nbdry_loc);
printf("\n\ncoords_loc: Lokale Koordinaten (crosspoints (immer 4), edge nodes, interior nodes):\n");
double* coords_loc = metra[domain]->domcoord;
for (size_t i=0;i<metra[domain]->ncoord_loc;i++){
   printf ("    (%lg,  %lg)\n", coords_loc[2*i], coords_loc[2*i+1] );
}
printf("\ndomelem: Elemente mit lokaler Knotennummerierung:\n");
for(size_t j=0; j<metra[domain]->nelem_loc; j++){
  printf("[");
  for(size_t k=0; k<3; k++){
    printf(" %td",metra[domain]->domelem[7*j+k]);
  }
  printf(" ]\n");
}
printf("\nc: Permutationsmatrix als Vektor:\n");
for(size_t i=0;i<metra[domain]->ncoord_loc;i++){
  printf("%td\n",metra[domain]->c[i]);
}
printf("\nnedgenodes (# Crosspoint+Edge Knoten = ) = %td", metra[domain]->nedgenodes);
printf("\nnfixed_loc (# Fixed Gebietsrandpunkte in Domain) = %td", metra[domain]->nfixed_loc);
printf("\n\nfixed_loc: Lokale fixed Knotennummern auf Gebietsrand:\n");
for(size_t i=0;i<metra[domain]->nfixed_loc;i++){
  printf("%td\n",metra[domain]->fixed_loc[i]);
}
printf("\nneighbours: Ranks of neighbours [s e n w] (-1 if none):\n");
for(size_t i=0;i<4;i++){
  printf("%td\n",metra[domain]->neighbours[i]);
}
printf("\nn_single_bdry: n_edgenodes (w.o. crosspoints) of each boundary sorted [s e n w]:\n");
for(size_t i=0;i<4;i++){
  printf("%td\n",metra[domain]->n_single_bdry[i]);
}
if(domain==0){
  printf("\nn_cross_glob (# parweise versch. crosspoints) = %td", metra[domain]->n_cross_glob);
  printf("\n\nc_cross: Ordnet allen crosspoints von allen ranks die globale Knotennummer zu:\n");
  for(size_t i=0;i<metra[domain]->n_cross_glob;i++){
    printf("%td\n",metra[domain]->c_cross[i]);
  }
}
printf("\nbool (1 = black, 0 = red) = %d", metra[domain]->black);

printf("\n------------------------\n\n");

return;
}
