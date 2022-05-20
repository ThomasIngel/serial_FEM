#ifndef _MESH_AUX_H
#define _MESH_AUX_H

#include "hpc.h"
#include "mesh_trans.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <blas_level1.h>

//Sortieralgorithmus
void sort_nodes(mesh_trans* metra,
		index *x, index n_coord, index n_coord_loc,
		index nedgenodes_sn, index nedgenodes_ew){
    index* sortet = metra -> c;
		       
    size_t count_interior = n_coord_loc-1;
    index count_south = 4;
    index count_east = 4 + nedgenodes_sn;
    index count_north = 4 + nedgenodes_sn + nedgenodes_ew;
    index count_west = 4 +  2 * nedgenodes_sn + nedgenodes_ew;

    for (size_t i = 0; i < n_coord; ++i){
	switch(x[i]) {
	    case 0: break;
	    case 1: break;
	    case 2:  sortet[count_interior--] = i; break;	//interior
	    case 21: sortet[count_south++] = i;  break;  	//south
	    case 16: sortet[count_east++] = i;  break;  	//east
	    case 25: sortet[count_north++] = i;  break;  	//north
	    case 11: sortet[count_west++] = i;  break;  	//west
	    //Crosspoints sortet (SE SW NW NE)
	    case 30: sortet[0] = i;  break;
	    case 35: sortet[1] = i;  break;
	    case 39: sortet[2] = i;  break;
	    case 34: sortet[3] = i;  break;
	}
    }	
}

//Initialisiert den Vektor c so, dass die Nummerierung der Knoten mit der,
//auf dem Aufgabenblatt uebereinstimmt,
//d. h. x = [x_Crosspoints; x_edges; x_interior]
//Dies ist notwendig fuer das parallele gs-verfahren
void sort_global_nodes(double* coord, index ncoord, index* c){

    //erste vier Knoten sind IMMER die Crosspoint nodes
    for(size_t i=0; i<4; i++){
	c[i] = i;
    }

    index count1 = 4;
    index count2 = ncoord-1;
    for(size_t i=4; i<ncoord; i++){
	if(coord[2*i] == 0 || coord[2*i] == 1 || 
	   coord[2*i+1] == 0 || coord[2*i+1] == 1){
	    c[count1] = i;
	    count1 += 1;
	}else{
	    c[count2] = i;
	    count2 -= 1;
	}
    }
}


//Gibt wieder, welche Elemente in der Domain mit den Knoten c sind
void get_dom_elements(index* elem, index nelem, index* c, index ncoord_loc,
	index nelem_loc, index* res){

    index gefunden = 1;
    index count;
    index count2 = 0;

    for(size_t i=0; i<nelem; i++){
	count = 0;
	gefunden = 1;
	for(size_t k=0; k<3; k++){
	    if (gefunden==1){
		for(size_t j=0; j<ncoord_loc; j++){
		    if(elem[7*i+k]==c[j]){
			count += 1;
			gefunden = 1;
			break;
		    }
		gefunden = 0;
		}
	    }
	}
	if (count==3){ //Element ist in Domain
	    res[count2] = i;
	    count2 += 1;
	}
	if(count2==nelem_loc) break;
    }
}


//Sucht v1 und v2 in vec
bool nodes_on_bdry(index v1, index v2, index* vec, index len){

    bool found1	=0;
    bool found2	=0;
    bool found	=0;

    for(size_t i=0;i<len;++i){
	if(vec[i] == v1) found1=1;
	if(vec[i] == v2) found2=1;
	if(found1 && found2){
	    found=1;
	    break;
	}
    }

    return found;
}


//Sucht target in einem Vektor c und gibt den Index zurueck
index search_index(index target, index* c, index len){

    for(index i=0;i<len;++i){
	if(c[i]==target){
	    return i;
	}
    }
    printf("search_index failed -> value %i not found\n",target);
}

//Suche Boundary Knoten/Arten fuer Domain
index* mesh_trans_getbdry(index* bdry, index nbdry,	//globale Daten aus mesh
	index* c, index ncoord_loc, index* nbdry_loc){	//lok. Daten mesh_trans

    index* dombdry = (index*) malloc(4*nbdry*sizeof(index));
    index anz = 0;

    for(size_t i=0; i<nbdry; i++){
	if(nodes_on_bdry(bdry[4*i], bdry[4*i+1], c, ncoord_loc)){
	    //Kante auf Bdr
	    for(size_t j=0; j<2; j++){
		dombdry[anz*4+j] = search_index(bdry[4*i+j], 
				    c, ncoord_loc);
		//Knoten mit lok. Knotennummer speichern
	    }
	    for(size_t j=2; j<4; j++){
		dombdry[anz*4+j] = bdry[4*i+j];
	    }
	    anz += 1;
	}
    }
    dombdry = realloc(dombdry, 4*anz*sizeof(index));
    *nbdry_loc = anz;

    return(dombdry);
}

//Suche die fixed Nodes in Teilgebiet
index* mesh_trans_getFixed(index* fixed, index nfixed,	//globale Daten aus mesh
	index* c, index ncoord_loc, index* nfixed_loc){	//lok. Daten mesh_trans

    index* fixed_loc = (index*) malloc(ncoord_loc*sizeof(index));
    index anz = 0;

    for(size_t i=0; i<ncoord_loc; i++){
	for(size_t j=0; j<nfixed; j++){
	    if(c[i]==fixed[j]){	    //Knoten ist auf Rand
		fixed_loc[anz] = i;
		anz += 1;
		break;
	    }
	}
    }
    fixed_loc = realloc(fixed_loc, anz*sizeof(index));
    *nfixed_loc = anz;

    return(fixed_loc);
}

//Findet die benachbarten Netze und setzt red black ordering
void 
get_neighbour_meshes(mesh_trans** metra, index i,index j, 
					    index m,index n){
    //south
    metra[i*n+j]->neighbours[0] = j==0   ? -1 : i*n + j - 1;
    //east
    metra[i*n+j]->neighbours[1] = i==m-1 ? -1 : (i+1)*n + j;
    //north
    metra[i*n+j]->neighbours[2] = j==n-1 ? -1 : i*n + j + 1;
    //west
    metra[i*n+j]->neighbours[3] = i==0   ? -1 : (i-1)*n + j;
    metra[i*n +j] -> black = (i%2) ^ (j%2) ==1 ? false : true;

}
void
prepare_crosspt_comm(mesh_trans** metra, index N){
    index* c_cross = (index*) malloc(4*N*sizeof(index));
    index* c_ = (index*) malloc(4*N*sizeof(index));
    index n_cross = 0;
    bool found = false;
    for (index i=0; i<N; ++i){
    	blasl1_icopy(metra[i]->c, &c_cross[i*4], 4, 1);
    }
    for (index i=0; i<4*N; ++i){
	for (index j =0; j<i; ++j){
	    if (c_cross[i] == c_cross[j]){
		found = true;
		c_[i] = c_[j];
		break;
	    }
	}
	if (!found){
	    c_[i] = n_cross++;
	}
    found = false;
    }
    metra[0] -> c_cross = c_cross;
    metra[0] -> cross_to_buf = c_;
    metra[0] -> n_cross_glob = n_cross;
}
//Diese Funktion splittet das ganze Netz in Teilnetze auf
void meshsplit(mesh* H, mesh_trans** metra, index anz_dom){

    index anz_dom_edgex, anz_dom_edgey;	//Anzahl der Domains in x-/y-Richtung

    //Gebietszerlegung in 2,4,8,16 domains moeglich
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

	printf("\nanz_dom_edgex = %td", anz_dom_edgex);
	printf("\nanz_dom_edgey = %td", anz_dom_edgey);
    index ncoord = H->ncoord;
	printf("\nncoord = %td", ncoord);
    double* coord = H->coord;
    index nelem  = H->nelem;
	printf("\nnelem = %td", nelem);
    index anz_nod_edge = sqrt(ncoord);											//Anzahl der nods in x/y-Richtung (Sqrt, da quadratisches Problem)
	printf("\nanz_nod_edge = %td", anz_nod_edge);
    index anz_nod_domedgex = (anz_nod_edge+anz_dom_edgex-1)/anz_dom_edgex;		//Anzahl nods pro domain in x
	printf("\nanz_nod_domedgex = %td", anz_nod_domedgex);
    index anz_nod_domedgey = (anz_nod_edge+anz_dom_edgey-1)/anz_dom_edgey;		//Anzahl nods pro domain in y
	printf("\nanz_nod_domedgey = %td", anz_nod_domedgey);
    index anz_nod_dom = anz_nod_domedgex*anz_nod_domedgey;						//Anzahl nods pro domain
	printf("\nanz_nod_dom = %td\n", anz_nod_dom);
    //Vektor soll wiedergeben ob Knoten sich x Wert in Domain befindet
    index* hilfsvek  = (index*) calloc(ncoord, sizeof(index));
    index* hilfsvek2  = (index*) calloc(ncoord, sizeof(index));

    //Vektr soll wiedergeben, welche Elemente in Domain sind
    index* elemindom = (index*) calloc((nelem/anz_dom),sizeof(index));

    //FINDE KNOTEN DIE ZUR DOMAIN GEHOEREN  (Domain Nummerierung z.B. 2 5 8
    //								      1 4 7
    //								      0 3 6

    //add values to determine position of boundarys and and crosspoint
    //x Koordinate
    for(size_t i=0; i<anz_dom_edgex; i++){
	for(size_t k=0; k<ncoord; k++){
	    if(coord[2*k]>=(double)i/anz_dom_edgex 
		    && coord[2*k]<= (double)(i+1)/anz_dom_edgex){						// Liegt x coord von Punkt in i-ter domain?
		if(coord[2*k]==(double)i/anz_dom_edgex){
		    hilfsvek[k] = 10; 	//west											// Linker Rand der domain
		} else if (coord[2*k]==(double)(i+1)/anz_dom_edgex){
		    hilfsvek[k] = 15;	//east											// Rechter Rand der domain
		} else {
		    hilfsvek[k] = 1;	//in einer der Domains							// Nicht am Rand --> in domain
		}
	    }else{
		hilfsvek[k] = 0;	//nicht in Domain
	    }
	}
	
	// -----> hilfsvek ist Länge ncoord und sagt für jeden nod, wo er sich in x-Richtung in seiner domain i befindet (Rand oder dazwischen)
	//y Koordinate
	for(size_t j=0; j<anz_dom_edgey; j++){
	    for(size_t k=0; k<ncoord; k++){
		hilfsvek2[k] = hilfsvek[k]; //copy vek (safe one loop)
		if(hilfsvek2[k]>0){ //Nur Pkte die in x Richtung in Frage kommen
		    if(coord[2*k+1]>=(double)j/anz_dom_edgey 
				&& coord[2*k+1]<=(double)(j+1)/anz_dom_edgey){
			if(coord[2*k+1]==(double)j/anz_dom_edgey){
			    hilfsvek2[k] += 20; 	//south
			} else if(coord[2*k+1]==(double)(j+1)/anz_dom_edgey){
			    hilfsvek2[k] += 24; 	//north
			} else{
			    hilfsvek2[k] += 1;	//in einer der Domains
			}
		    }
		}else{
		}
	    }
	    metra[i*anz_dom_edgey+j]->nedgenodes =  2*(anz_nod_domedgex-1)+
						    2*(anz_nod_domedgey-1);
	    metra[i*anz_dom_edgey+j]->ncoord_loc = anz_nod_dom;
	    metra[i*anz_dom_edgey+j]->ncoord_glo = ncoord;
		metra[i*anz_dom_edgey+j]->nelem_loc  = nelem/anz_dom;
	    sort_nodes(metra[i*anz_dom_edgey+j],hilfsvek2, ncoord, 
		       anz_nod_dom,anz_nod_domedgex-2, anz_nod_domedgey-2);
	    
	    metra[i*anz_dom_edgey+j]->n_single_bdry[0] = anz_nod_domedgex-2;
	    metra[i*anz_dom_edgey+j]->n_single_bdry[1] = anz_nod_domedgey-2;
	    metra[i*anz_dom_edgey+j]->n_single_bdry[2] = anz_nod_domedgex-2;
	    metra[i*anz_dom_edgey+j]->n_single_bdry[3] = anz_nod_domedgey-2;

	    get_neighbour_meshes(metra, i, j,anz_dom_edgex, anz_dom_edgey);
	}
    }

    //UEBERNEHME COORDINATES UND ELEMENTS MIT NEUER LOKALER KNOTENNUMMERIERUNG
    for(size_t i=0; i<anz_dom; i++){
	index* c = metra[i]->c;
	for(size_t j=0; j<anz_nod_dom; j++){
	    metra[i]-> domcoord[2*j] = coord[2*c[j]];	    //Koordinates anleg.
	    metra[i]-> domcoord[2*j+1] = coord[2*c[j]+1];
	}
	//Coordinates mit neuen Knotennummerierung
	get_dom_elements(H->elem, nelem, metra[i]->c, metra[i]->ncoord_loc,
		metra[i]->nelem_loc, elemindom);
	
	for(size_t j=0; j<metra[i]->nelem_loc; j++){
	    for(size_t k=0; k<3; k++){
		metra[i]->domelem[7*j+k] =
		    search_index(H->elem[7*elemindom[j]+k], metra[i]->c,
			    metra[i]->ncoord_loc);
	    }
	}
    }

    //get bdry and fixed nodes for each rank
    for(size_t i=0; i<anz_dom; i++){
	metra[i]->dombdry =  mesh_trans_getbdry(H->bdry, H->nbdry,
		    metra[i]->c, metra[i]->ncoord_loc, &metra[i]->nbdry_loc);
	metra[i]->fixed_loc =  mesh_trans_getFixed(H->fixed, H->nfixed,
		    metra[i]->c, metra[i]->ncoord_loc, &metra[i]->nfixed_loc);
    }

    //prepare communication
    prepare_crosspt_comm(metra,anz_dom);

    free(hilfsvek);
    free(hilfsvek2);
    free(elemindom);
}

#endif

