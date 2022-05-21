/**
 * @file sed_sm_build.c
 * @brief Funktionen zum aufstellen der Steifigkeitsmatrix.
 * @author Joshua Bauske, Heiko Karus, Ernst Tobias Leutbecher, Benjamin Fechtig, Felix
 * Götz
 * @date 17.09.2021
 */

#include "hpc.h"
#include "mesh_trans.h"

/**
 * Funktion zur Berechnung des SED Matrix Speicherlayouts.
 *
 * @param[in] mesh_loc	Gitter des lokalen Teilgebiets als mesh_trans struct.
 *
 * @return Zurückgeben einer allokierten SED Matrix.
 */
sed *sed_sm_pattern(mesh_trans *mesh_loc)
{
	// Verschiedene Variablen und Zeiger für die Berechnungen
	index k, j, n, p, nC, nT, nE, nz, *Elem, ind[3], *Si, *w, imin, imax;
	sed *S;

	// Indizes in x und y Richtung wie die Steifigkeitsmatrix für ein Element ausschaut
	static int ai[3] = {0, 0, 1}, aj[9] = {1, 2, 2};

	// Auslesen der Daten vom mesh_trans Objekt
	nT = mesh_loc->nelem_loc;
	nC = mesh_loc->ncoord_loc;
	Elem = mesh_loc->domelem;

	// get structure of sparse matrix
	n = nC; // Dimension der Steifigkeitsmatrix
	nz = n + 1 + 3 * nT;

	// Speicher anlegen für die Matrix S und den temporären Speicher w
	S = sed_alloc(n, nz, 0);
	Si = S->i;
	w = calloc(n, sizeof(index));
	if (!S || !w) {
		return sed_done(S, w, NULL, 0); // out of memory
	}

	// column counts
	for (k = 0; k < nT; k++) {
		for (j = 0; j < 3; j++) {
			ind[j] = Elem[7 * k + j]; // Knotennummer des Elements k erhalten
		}
		for (j = 0; j < 3; j++) {
			// Speichern, wie viele Nebendiagonaleinträge es in jeder Spalte gibt
			// Symmetrie wird direkt mit berücksichtigt!!!
			w[HPC_MIN(ind[ai[j]], ind[aj[j]])]++;
		}
	}

	// Berechnen der Zeiger für die Spalten:
	//  - Erst Kummulative Summe bilden,
	//  - dann die Anzahl n Diagonalelementen für den Spalten Pointer aufaddieren
	hpc_cumsum(Si, w, n);
	for (k = 0; k < n; k++) {
		w[k] += n + 1;
	}
	for (k = 0; k < n + 1; k++) { // An n kommt die Anzahl an Einträgen
		Si[k] += n + 1;
	}
	// Berechnen der Zeilen indizes:
	for (k = 0; k < nT; k++) {
		// Knotennummern des Elements k bestimmen
		for (j = 0; j < 3; j++) {
			ind[j] = Elem[7 * k + j];
		}
		// Bestimmen des Zeilen indizes vom entprechenden Knoten und Speichern im Indize
		// Zeiger Si
		for (j = 0; j < 3; j++) {
			imin = HPC_MIN(ind[ai[j]], ind[aj[j]]);
			imax = HPC_MAX(ind[ai[j]], ind[aj[j]]);
			Si[w[imin]++] = imax;
		}
	}

	// Freigeben des Workpace Speichers und löschen doppelter Einträge
	free(w);
	if (!sed_dupl(S)) {
		return sed_free(S);
	}

	return S;
}

/**
 * Funktion zur Berechnung der Steifigkeitsmatrix eines Elements.
 * ordering w.r.t. [ p1, p2, p3, m1=(p1+p2)/2, m2=(p2+p3)/2, m3=(p3+p1)/2]
 *
 * @param[in] p1 Koordinaten Knoten 1.
 * @param[in] p2 Koordinaten Knoten 2.
 * @param[in] p3 Koordinaten Knoten 3.
 * @param[in] dx Speicher für die Diagonalelemente.
 * @param[in] ax Speicher für die Nebendiagonalelemente.
 */
void sed_sm_element(double p1[2], double p2[2], double p3[2], double dx[3], double ax[3])
{
	// Flächeninhalt des gesamten Elements (triangle) T berechnen
	double T =
		0.5 * ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]));

	// Faktor für Berechnungsformel der stiffness Matrix Einträge
	double fac = 1 / (4 * T);

	// Mehrfach verwendete zwischenergebnis
	double delta[4];
	delta[0] = (p2[1] - p3[1]); // y2-y3
	delta[1] = (p3[0] - p2[0]); // x3-x2
	delta[2] = (p3[1] - p1[1]); // y3-y1
	delta[3] = (p1[0] - p3[0]); // x1-x3

	// Steifigkeitsmatrix m als dx (Diagonale) and ax (Nebendiagonale)
	dx[0] = fac * (delta[0] * delta[0] + delta[1] * delta[1]); // S_11
	ax[0] = fac * (delta[0] * delta[2] + delta[1] * delta[3]); // S_21 = S_12
	ax[1] = -dx[0] - ax[0];									   // S_31 = S_13
	dx[1] = fac * (delta[2] * delta[2] + delta[3] * delta[3]); // S_22
	ax[2] = -ax[0] - dx[1];									   // S_32 = S_23
	dx[2] = -ax[1] - ax[2];									   // S_33
}

/**
 * Funktion zur Berechnung des SED Matrix Speicherlayouts.
 *
 * @param mesh_loc  Lokales Gitter als mesh_trans.
 *
 * @return Zurückgeben einer allokierten gefüllten SED Matrix.
 */
sed *sed_sm_build(mesh_trans *mesh_loc)
{
	// Verschiedene Variablen und Zeiger für die Berechnungen
	index j, k, n, p, nC, nT, nz, *Elem, ind[3], *Ai, *w, imin, imax;
	double dx[3], ax[3], *Coord, *Ax;
	sed *A;

	// Indizes in x und y Richtung wie die Steifigkeitsmatrix für ein Element
	// ausschaut
	static int ai[3] = {0, 0, 1}, aj[9] = {1, 2, 2};

	// Auslesen der Daten vom mesh Objekt
	nT = mesh_loc->nelem_loc;
	nC = mesh_loc->ncoord_loc;
	Coord = mesh_loc->domcoord;
	Elem = mesh_loc->domelem;

	// Auslesen der Daten vom sed Objekt und Speicher anlegen für die Matrixwerte
	A = sed_sm_pattern(mesh_loc);
	n = A->n;
	Ai = A->i;
	if (!(A->x)) {
		Ax = A->x = calloc(Ai[n], sizeof(double)); // Ai[n] = A->nzmax
	}
	if (!Ax) {
		return (0);
	}

	// Für jedes Element die Steifigkeitsmatrix berechnen
	for (k = 0; k < nT; k++) {
		// Knotennummer des Elements k erhalten
		for (j = 0; j < 3; j++) {
			ind[j] = Elem[7 * k + j];
		}

		// Berechnen der Steifigkeitsmatrix des Elements k
		sed_sm_element(Coord + 2 * ind[0], // Koordinaten Knoten 1
					   Coord + 2 * ind[1], // Koordinaten Knoten 2
					   Coord + 2 * ind[2], // Koordinaten Knoten 3
					   dx,				   // Diagonalelemente
					   ax				   // Nebendiagonalelemente
		);

		// Füllen der gesammten SED Steifigkeitsmatrix
		for (j = 0; j < 3; j++) {
			Ax[ind[j]] += dx[j]; // Einsetzen der Diagonalelemente
		}
		for (j = 0; j < 3; j++) {
			imin = HPC_MIN(ind[ai[j]], ind[aj[j]]);
			imax = HPC_MAX(ind[ai[j]], ind[aj[j]]);
			for (p = Ai[imin]; p < Ai[imin + 1]; p++) {
				if (Ai[p] == imax) {
					Ax[p] += ax[j]; // Einsetzen der Nichtdiagonal Elemente
					break;
				}
			}
		}
	}

	return A;
}