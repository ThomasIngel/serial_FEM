#include "hpc.h"
#include "mesh_trans.h"

/**
 * Funktion zur Berechnung der rechten Seite der Neumann-Randbedingung.
 *
 * @param[in] p1    Koordinaten Knoten 1.
 * @param[in] p2    Koordinaten Knoten 2.
 * @param[in] typ   Gebietsnummer, falls unterschiedliche Gebiete vorhanden.
 * @param[in] fc    Funktion für die Neumann-Randbedingung.
 * @param[out] b    Berechnete Funktionswerte für Knoten 1 und 2.
 */
void rhs_Neumann(double p1[2], double p2[2], index typ, double (*fc)(double *, index),
				 double b[2])
{
	// Verschiedene Variablen und Zeiger für die Berechnungen
	int i;
	double x[2], h = 0.0, fac;
	// Mittelpunkts-Regel
	for (i = 0; i < 2; i++) {
		// Euklidische Norm berechnen 1.
		h += (p2[i] - p1[i]) * (p2[i] - p1[i]);
		// Koordinaten Mittelpunkt zwischen Startknoten und Endknoten
		x[i] = (p1[i] + p2[i]) / 2.0;
	}
	// Euklidische Norm berechnen 2. = Distanz zwischen Startknoten und Endknoten
	h = sqrt(h);
	// Distanz bis zur Mitte benötigt: * 0.5
	h = h / 2.0;

	fac = fc(x, typ) * h;
	b[0] = fac;
	b[1] = fac;
}

/**
 * Funktion zur Berechnung der rechten Seite der Volumenkraft.
 *
 * @param[in] p1    Koordinaten Knoten 1.
 * @param[in] p2    Koordinaten Knoten 2.
 * @param[in] typ   Randbedingungs Nummmer
 * @param[in] fc    Funktion für die Neumann-Randbedingung.
 * @param[out] b    Berechnete Funktionswerte für Knoten 1 und 2.
 */
void rhs_Volumen(double p1[2], double p2[2], double p3[2], index typ,
				 double (*fc)(double *, index), double b[3])
{
	// Verschiedene Variablen und Zeiger für die Berechnungen
	int i;
	double mid[2], d[2][2], fac;

	// Mittelpunkts-Regel
	for (i = 0; i < 2; i++) {
		d[0][i] = p1[i] - p3[i];
		d[1][i] = p2[i] - p1[i];
		mid[i] = (p1[i] + p2[i] + p3[i]) / 3.0;
	}
	fac = fc(mid, typ) / 6.0 * (d[0][0] * d[1][1] - d[1][0] * d[0][1]);
	for (i = 0; i < 3; i++) {
		b[i] = fac;
	}
}

/**
 * Funktion zur Berechnung der rechten Seite
 *
 * @param[in] mesh_loc  Lokales Teilgebiet als mesh_transfer Objekt.
 * @param[in] b         Vektor der rechten Seite des LGS.
 * @param[in] fV        Funktion der Volumenkraft f.
 * @param[in] fN        Funktion der Neumann-Randbedingung.
 */
void mesh_trans_rhs(const mesh_trans *mesh_loc, double *b,
					double (*fV)(double *, index), double (*fN)(double *, index))
{
	// Verschiedene Variablen und Zeiger für die Berechnungen
	index j, k, nT, nB, *Elem_loc, *Bdry_loc, ind[3];
	double *Coord_loc, bx[3];

	// Auslesen der Daten vom mesh Objekt
	nT = mesh_loc->nelem_loc;
	nB = mesh_loc->nbdry_loc;
	Coord_loc = mesh_loc->domcoord;
	Elem_loc = mesh_loc->domelem;
	Bdry_loc = mesh_loc->dombdry;

	// Rechte Seite: Volumenkraft f
	for (k = 0; k < nT; k++) {
		// Knotennummer des Elements k erhalten
		for (j = 0; j < 3; j++) {
			ind[j] = Elem_loc[7 * k + j];
		}

		rhs_Volumen(Coord_loc + 2 * ind[0], // Koordinaten Knoten 1
					Coord_loc + 2 * ind[1], // Koordinaten Knoten 2
					Coord_loc + 2 * ind[2], // Koordinaten Knoten 3
					Elem_loc[7 * k + 6],	// Subdomain nummer
					fV,						// Volumenkraft Funktion
					bx						// Funktionswerte der Knoten
		);

		// 3 Knoten pro Element
		for (j = 0; j < 3; j++) {
			b[ind[j]] += bx[j];
		}
	}

	// Rechte Seite: Neumann-Randbedingung
	for (k = 0; k < nB; k++) {
		if (Bdry_loc[4 * k + 3]) {
			// Knotennummern der Kante
			for (j = 0; j < 2; j++) {
				ind[j] = Bdry_loc[4 * k + j];
			}

			rhs_Neumann(Coord_loc + 2 * ind[0], // Koordinaten Knoten 1
						Coord_loc + 2 * ind[1], // Koordinaten Knoten 2
						Bdry_loc[4 * k + 3],	// Randbedingungs Nummmer
						fN,						// Neumanrand Funktion
						bx						// Funktionswerte der Knoten
			);

			// 2 Knoten pro Kante
			for (j = 0; j < 2; j++) {
				b[ind[j]] += bx[j];
			}
		}
	}
}
