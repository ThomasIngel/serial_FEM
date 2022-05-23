#ifndef _SCATTER_MESH_
#define _SCATTER_MESH_

#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <stdbool.h>
#include "hpc.h"
#include "mesh_trans.h"

//Quelle: https://stackoverflow.com/questions/40807833/sending-size-t-type-data-with-mpi
#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif

#define numuni 3    //Anzahl der auf allen Ranks gleich langen Vektoren
#define numnotuni 2 //Anzahl der nicht gleich langen Vektoren
#define single 9    //Anzahl an einzelnen Werten

/*//---------------------------------KLASSE---------------------------------------
typedef struct mesh_transfer_data{
    1 index ncoord_loc;		//Anzahl Koordinaten lokal
    2 index ncoord_glo;		//Anzahl Koordinaten global
    3 index nelem_loc;		//Anzahl Elemente lokal
    4 index nbdry_loc;		//Anzahl Gebietsrandpunkte im Gebiet

    9 double* domcoord;		//lokale Koordinaten
    7 index* domelem;		//Elemente, mit lokaler Knotennummerierung
    8 index* c;			//C Matrix als Vektor
    10 index* dombdry;		//Knoten am Rand und Art der Randbed.

    5 index nedgenodes;		//Anzahl Crosspoint+Edge Knoten -> Es gibt immer
				//vier Crosspoints! Der domcoord Vektor ist ge-
				//ordnet nach
				//[crosspoints   edge nodes   interior nodes]'

    6 index nfixed_loc;		//Anzahl der Gebietsrandpunkte in Domain
    11 index* fixed_loc;	//Lokale Knotennummern der Knoten, welche auf
				//dem Gebietsrand liegen
    12 index neighbours[4]; 	//ranks of the neighbours [s e n w] if no neigh-
    				//bour available -1
    13 index n_single_bdry[4];  //n_edgenodes (w.o. crosspoints) of each boundary 
    				//sortet [s e n w]
    14 bool black;			//meshes have a red black ordering, true is 
    				//black, false is red
} mesh_trans ;*/

mesh_trans* scatter_meshes(mesh_trans** global_mesh,MPI_Comm comm,
	index domains, index dof){

	int rank; int size;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&size);
	int worker=size-1;

	if(rank==0){

	    //Zuerst noch fixed_loc und dombdry allokieren
	    int req_len=worker*(single+numuni+numnotuni);
	    MPI_Request request[req_len];
	    int count=0;
	    //Sende parallel mit Isend zunaechst die einzelnen Werte
	    for(size_t i=1; i<size; ++i){
					
		//Sende ncoord_loc
		MPI_Isend(&global_mesh[i]->ncoord_loc,1,my_MPI_SIZE_T,i,1,
			comm,&request[count++]);

		//Sende ncoord_glo
		MPI_Isend(&global_mesh[i]->ncoord_glo,1,my_MPI_SIZE_T,i,2,
			comm,&request[count++]);

		//Sende nelem_loc
		MPI_Isend(&global_mesh[i]->nelem_loc,1,my_MPI_SIZE_T,i,3,
			comm,&request[count++]);

		//Sende nbry_loc 
		MPI_Isend(&global_mesh[i]->nbdry_loc,1,my_MPI_SIZE_T,i,4,
			comm,&request[count++]);
		
		//Sende nedgenodes
		MPI_Isend(&global_mesh[i]->nedgenodes,1,my_MPI_SIZE_T,i,5,
			comm,&request[count++]);

		//Sende nfixed_loc
		MPI_Isend(&global_mesh[i]->nfixed_loc,1,my_MPI_SIZE_T,i,6,
			comm,&request[count++]);
		
		//Sende neighbours
		MPI_Isend(global_mesh[i]->neighbours,4,my_MPI_SIZE_T,i,12,
		          comm, &request[count++]);
		//Sende n_single_bdry
		MPI_Isend(global_mesh[i]->n_single_bdry,4,my_MPI_SIZE_T,i,13,
		          comm, &request[count++]);
		//Sende black
		MPI_Isend(&global_mesh[i]->black,1,MPI_C_BOOL,i,14,
			  comm, &request[count++]);
		//Sende domelem Vektor 
		MPI_Isend(global_mesh[i]->domelem, 7*global_mesh[i]->nelem_loc,
			my_MPI_SIZE_T,i,7,comm,&request[count++]);

		//Sende c Vektor 
		MPI_Isend(global_mesh[i]->c, global_mesh[i]->ncoord_loc,
			my_MPI_SIZE_T,i,8,comm,&request[count++]);
		
		//Sende domcoord Vektor 
		MPI_Isend(global_mesh[i]->domcoord,2*global_mesh[i]->ncoord_loc,
			MPI_DOUBLE,i,9,comm,&request[count++]);

		//Versende Daten der non_uniform Pointer 
		MPI_Isend(global_mesh[i]->dombdry,4*global_mesh[i]->nbdry_loc,
			    my_MPI_SIZE_T,i,10,comm,&request[count++]);
		MPI_Isend(global_mesh[i]->fixed_loc,global_mesh[i]->nfixed_loc,
			    my_MPI_SIZE_T,i,11,comm,&request[count++]);

	    }
	    //Synchronisiere, zwingend notwendig bevor es weitergeht
	    for(size_t i=0;i<req_len;++i){
		MPI_Wait(&request[i],MPI_STATUS_IGNORE);
	    }
	    return global_mesh[0];
	}else{
	    mesh_trans* local_mesh=alloc_mesh_trans(domains,dof);
	    MPI_Request request_s[single];
	    int count=0;

	    //Normale Skalare
	    MPI_Irecv(&local_mesh->ncoord_loc,1,my_MPI_SIZE_T,0,1,
		    comm,&request_s[count++]);

	    MPI_Irecv(&local_mesh->ncoord_glo,1,my_MPI_SIZE_T,0,2,
		    comm,&request_s[count++]);

	    MPI_Irecv(&local_mesh->nelem_loc,1,my_MPI_SIZE_T,0,3,
		    comm,&request_s[count++]);

	    MPI_Irecv(&local_mesh->nbdry_loc,1,my_MPI_SIZE_T,0,4,
		    comm,&request_s[count++]);

	    MPI_Irecv(&local_mesh->nedgenodes,1,my_MPI_SIZE_T,0,5,
		    comm,&request_s[count++]);

	    MPI_Irecv(&local_mesh->nfixed_loc,1,my_MPI_SIZE_T,0,6,
		    comm,&request_s[count++]);

	    MPI_Irecv(local_mesh->neighbours,4,my_MPI_SIZE_T,0,12,
		    comm,&request_s[count++]);

	    MPI_Irecv(local_mesh->n_single_bdry,4,my_MPI_SIZE_T,0,13,
		    comm,&request_s[count++]);

	    MPI_Irecv(&local_mesh->black,1,MPI_C_BOOL,0,14,
		    comm,&request_s[count++]);
	    for(size_t i=0;i<single;++i){
		MPI_Wait(&request_s[i],MPI_STATUS_IGNORE);
	    }

	    MPI_Request request_u[numuni];
	    count=0;
	    //Uniform Pointer
	    MPI_Irecv(local_mesh->domelem,7*local_mesh->nelem_loc,
		    my_MPI_SIZE_T,0,7, comm,&request_u[count++]);

	    MPI_Irecv(local_mesh->c,local_mesh->ncoord_loc,my_MPI_SIZE_T,0,8,
		    comm,&request_u[count++]);

	    MPI_Irecv(local_mesh->domcoord,2*local_mesh->ncoord_loc,
		    MPI_DOUBLE,0,9, comm,&request_u[count++]);

	    for(size_t i=0;i<numuni;++i){
		MPI_Wait(&request_u[i],MPI_STATUS_IGNORE);
	    }
	    //lege die Speicherflaechen auf dem Heap an
	    local_mesh->dombdry=malloc(sizeof(index)*4*local_mesh->nbdry_loc);
	    local_mesh->fixed_loc=
		malloc(sizeof(index)*local_mesh->nfixed_loc);

	    MPI_Request request_nu[numnotuni];
	    count=0;
	    //Non uniform pointers empfangen
	    MPI_Irecv(local_mesh->dombdry,4*local_mesh->nbdry_loc,
		    my_MPI_SIZE_T,0,10, comm,&request_nu[count++]);

	    MPI_Irecv(local_mesh->fixed_loc,local_mesh->nfixed_loc,
		    my_MPI_SIZE_T,0,11, comm,&request_nu[count++]);

	    //Synchronisiere
	    for(size_t i=0;i<numnotuni;++i){
		MPI_Wait(&request_nu[i],MPI_STATUS_IGNORE);
	    }
	    return local_mesh;
	}
	
}

#endif





