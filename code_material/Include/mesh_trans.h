#ifndef _MESH_TRANS_H
#define _MESH_TRANS_H

#include <mpi.h>
#include "hpc.h"
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

/**
 * @brief Kommunikationsklasse, mit der das Mesh upgedated wird.
 * @file mesh_trans.h
 * 
 */

//---------------------------------KLASSE---------------------------------------
/**
 * @brief Enthält Daten, die für die Kommunikation zwischen den einzelnen Domains
 *        notwendig ist, sowie das lokale mesh.
 * 
 */
typedef struct mesh_transfer_data {
    index ncoord_loc;                /** @brief Anzahl Koordinaten lokal */
    index ncoord_glo;                /** @brief Anzahl Koordinaten global */
    index nelem_loc;                 /** @brief Anzahl Elemente lokal */
    index nbdry_loc;                 /** @brief Anzahl Gebietsrandpunkte im Gebiet */

    double* domcoord;                /** @brief lokale Koordinaten */
    index* domelem;                  /** @brief Elemente, mit lokaler Knotennummerierung */
    index* c;                        /** @brief C Permutationsmatrix als Vektor */
    index* dombdry;                  /** @brief Knoten am Rand und Art der Randbed. */

    index nedgenodes;                /** @brief Anzahl Crosspoint+Edge Knoten -> Es gibt immer
                                         vier Crosspoints! Der domcoord Vektor ist ge-
                                         ordnet nach [crosspoints   edge nodes   interior nodes]' */

    index nfixed_loc;                /** @brief Anzahl der Gebietsrandpunkte in Domain */
    index* fixed_loc;                /** @brief Lokale Knotennummern der Knoten, welche auf
                                         dem Gebietsrand liegen */

    
    index neighbours[4];    /** @brief ranks of the neighbours [s e n w] if no neigh-
                                bour available -1 */
    index n_single_bdry[4]; /** @brief n_edgenodes (w.o. crosspoints) of each boundary 
                                    sortet [s e n w] */
    index n_cross_glob;     /** @brief nur auf rank 0, anzahl paarweise vers crosspoints */
    index* c_cross;         /** @brief nur auf rank 0, ordnet allen crosspoints von a
                                allen ranks die globale Knotennummer zu */
    index* cross_to_buf;    /** @brief nur auf rank 0, index fuer den buffer zur 
                                addition der crosspoints auf rank 0 */
    double* buf_cross;      /** @brief nur auf rank 0 buffer fuer gather crosspoints
                                Laenge des buffers ist 4*size */
    double* buf_cross_aggr; /** @brief nur auf rank 0 buffer um crosspoints 
                                aufzuaddieren, Laenge ist n_crosspoints_global */
                                
    double* buf_bdry[4];    /** @brief buffer um werte der boundaries empfangen zu 
                                koennen, 2_dimensionales array buf_brdy[0] fuer 
                                south bdry [1] fuer east usw. */
    bool black;             /** @brief meshes have a red black ordering, true is 
                                black, false is red */
    /* probably unneccessary
    index* map_south;
    index* map_east;
    index* map_north;
    index* map_west;
    */
} mesh_trans ;

mesh_trans* alloc_mesh_trans(index anz_dom, index dof);
mesh_trans* free_mesh_trans(mesh_trans* metra);
void sort_nodes(mesh_trans* metra,index *x, index n_coord, index n_coord_loc,index nedgenodes_sn, index nedgenodes_ew);
void sort_global_nodes(double* coord, index ncoord, index* c);
void get_dom_elements(index* elem, index nelem, index* c, index ncoord_loc,index nelem_loc, index* res);
bool nodes_on_bdry(index v1, index v2, index* vec, index len);
index search_index(index target, index* c, index len);
index* mesh_trans_getbdry(index* bdry, index nbdry,index* c, index ncoord_loc, index* nbdry_loc);
index* mesh_trans_getFixed(index* fixed, index nfixed,index* c, index ncoord_loc, index* nfixed_loc);
void get_neighbour_meshes(mesh_trans** metra, index i,index j,index m,index n);
void prepare_crosspt_comm(mesh_trans** metra, index N);
void meshsplit(mesh* H, mesh_trans** metra, index anz_dom);
void mesh_trans_print (mesh_trans** metra, index domain);

// mesh_buildRhs_loc(const mesh_trans *M, double *b, double (*fV)(double *, index),double (*fN)(double *, index));

sed *sed_sm_pattern(mesh_trans *mesh_loc);
void sed_sm_element(double p1[2], double p2[2], double p3[2], double dx[3], double ax[3]);
sed *sed_sm_build(mesh_trans *mesh_loc);

mesh_trans* scatter_meshes(mesh_trans** global_mesh,MPI_Comm comm,index domains, index dof);

void mesh_trans_rhs(const mesh_trans *mesh_loc, double *b,double (*fV)(double *, index), double (*fN)(double *, index));
void rhs_Volumen(double p1[2], double p2[2], double p3[2], index typ,double (*fc)(double *, index), double b[3]);
void rhs_Neumann(double p1[2], double p2[2], index typ, double (*fc)(double *, index), double b[2]);




#endif
