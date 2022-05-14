#include "mesh_trans.h"
#include "hpc.h"

mesh_trans* alloc_mesh_trans(index anz_dom, index dof){
    index anz_dom_edgex, anz_dom_edgey;//Azahl der Domains in x-/y-Richtung
    //Gebietszerlegung nur in 2,4,8,16 domains möglich
    if(anz_dom == 2){
	anz_dom_edgex = 2;
	anz_dom_edgey = 1;
    }else if(anz_dom == 4){
	anz_dom_edgex = 2;
	anz_dom_edgey = 2;
    }else if(anz_dom == 8){
	anz_dom_edgex = 4;
	anz_dom_edgey = 2;
    }else if(anz_dom == 16){
	anz_dom_edgex = 4;
	anz_dom_edgey = 4;
    }else{
	printf("Keine Gebietszerlegung mit angegeber Domainanzahl moeglich\n");
	abort();
    }

    mesh_trans* metra = (mesh_trans*) malloc(sizeof(mesh_trans));
    index dof_edge = sqrt(dof);
    //index anz_dom_edge = sqrt(anz_dom);
    index nodes_per_sidex = (dof_edge+anz_dom_edgex-1)/anz_dom_edgex;
    index nodes_per_sidey = (dof_edge+anz_dom_edgey-1)/anz_dom_edgey;

    metra->domcoord = (double*) 
	malloc(2*nodes_per_sidex*nodes_per_sidey*sizeof(double));
    metra->domelem  = (index*) malloc(7*(nodes_per_sidex-1)*(nodes_per_sidey-1)*
					    2*sizeof(index));
    metra->c = (index*)malloc(nodes_per_sidex*nodes_per_sidey*sizeof(index));
    
    metra->fixed_loc = NULL;
    metra->dombdry =  NULL;
    metra->c_cross = NULL;
    metra->cross_to_buf = NULL;
    metra->buf_cross = NULL;
    metra->buf_cross_aggr = NULL;
    metra->buf_bdry[0]= NULL;
    metra->buf_bdry[1]= NULL;
    metra->buf_bdry[2]= NULL;
    metra->buf_bdry[3]= NULL;

    /*
    metra->map_south =  NULL;
    metra->map_east =  NULL;
    metra->map_north =  NULL;
    metra->map_west =  NULL;
    */
    return metra;
}

mesh_trans* free_mesh_trans(mesh_trans* metra){
    if (!metra) return (NULL);
    free(metra->domcoord);
    free(metra->domelem);
    free(metra->c);
    if(metra->fixed_loc) free(metra->fixed_loc);
    if(metra->dombdry) free(metra->dombdry);
    if(metra->c_cross) free(metra->c_cross);
    if(metra->cross_to_buf) free(metra->cross_to_buf);
    if(metra->buf_cross) free(metra->buf_cross);
    if(metra->buf_cross_aggr) free(metra->buf_cross_aggr);
    if(metra->buf_bdry[0]) free(metra->buf_bdry[0]);
    if(metra->buf_bdry[1]) free(metra->buf_bdry[1]);
    if(metra->buf_bdry[2]) free(metra->buf_bdry[2]);
    if(metra->buf_bdry[3]) free(metra->buf_bdry[3]);
    /*
    if(metra->map_south) free(metra->map_south);
    /*
    if(metra->map_south) free(metra->map_south);
    /*
    if(metra->map_south) free(metra->map_south);
    /*
    if(metra->map_south) free(metra->map_south);
    if(metra->map_east) free(metra->map_east);
    if(metra->map_north) free(metra->map_north);
    if(metra->map_west) free(metra->map_west);
    */
    free(metra);
    return(NULL);
}
