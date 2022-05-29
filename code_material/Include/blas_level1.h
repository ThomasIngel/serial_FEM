#ifndef BLAS_LEVEL1_H
#define BLAS_LEVEL1_H

#include <stddef.h>
#include <stdio.h>

/** 
 * @file blas_level1.h
 * @brief Blas Level1 Vektor-Funktionen
 * 
 * Syntax allgemein: blasl1_<typ><op>(....) 
 * e.g: blasl1_ddot() -> d = double, dot = skalar produkt
 * e.g: blasl1_iprint -> i= index, print = print vektor
 */


/** 
 * @brief Skaliert Vektor x<- alpha*x
 *
 * @param x double Vektor
 * @param length Laenge des Vektors
 * @param alpha double Skalierungsfaktor
 * @retrun überschreibt Vektor x
 */
void
blasl1_dscal(double * x, size_t len, const double alpha);


/** 
 * @brief Skalarprodukt res = x^T * y
 *
 * @param x double Vektor
 * @param y double Vektor
 * @param length Laenge der Vektoren
 * @retrun skalarprodukt der Vektoren
 */
double
blasl1_ddot(const double * x, const double * y, size_t len);

/**
 * @brief Vektor-Update x = beta * x + alpha * y
 * 
 * @param x Vektor
 * @param y Vektor
 * @param length Länge der Vektoren x und y
 * @param alpha 
 * @param beta 
 */
void
blasl1_daxpy(double* x, double* y, index length, double alpha, double beta);
// copy Functions

/**
 * @brief Double Vektor kopieren b = alpha * a
 * 
 * @param a Vektor
 * @param b Vektor
 * @param length Länge von a und b
 * @param alpha 
 */
void
blasl1_dcopy(const double* a, double* b, index length, double alpha);

/**
 * @brief Index Vektor kopieren b = alpha * a
 * 
 * @param a Vektor
 * @param b Vektor
 * @param length Länge von a und b
 * @param alpha 
 */
void
blasl1_icopy(const index* a, index* b, index length, double alpha);

void blasl1_intcopy(const int* a, int* b, index length, double alpha);

// print functions

/**
 * @brief Print Double-Vektor, jeder Wert in 1 Zeile
 * 
 * @param data Vektor
 * @param len Länge des Vektors
 */
void blasl1_dprint(double * data, size_t len);

/**
 * @brief Print index-Vektor, jeder Wert in 1 Zeile
 * 
 * @param data Vektor
 * @param len Länge des Vektors
 */
void blasl1_iprint(index* data, size_t len);

void blasl1_intprint(index* data, size_t len);


#endif //BLAS_LEVEL1_H
