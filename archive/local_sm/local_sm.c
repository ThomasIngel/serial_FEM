#include "hpc.h"
#include "mesh_trans.h"

int main() {

  mesh* H = get_refined_mesh(1);
  mesh_print(H,0);

  index anz_dom = 2;
  index ncoord = H->ncoord;
  mesh_trans **metra;
  metra = malloc ( (anz_dom) * sizeof(mesh_trans));

  for(size_t i=0;i<anz_dom;i++){
    metra[i]=alloc_mesh_trans(anz_dom,ncoord);
  }

  meshsplit(H, metra, anz_dom);
  
  /*index domain = 0;
  mesh_trans_print (metra, domain);*/

  sed** S;
  S = malloc ( (anz_dom) * sizeof(sed));
  for (size_t i=0;i<anz_dom;i++){
    S[i] = sed_sm_build(metra[i]);
    sed_print(S[i],0);
  }
  

}
