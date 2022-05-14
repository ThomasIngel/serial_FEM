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
blasl1_dscal(double * x, size_t len, const double alpha){
    if (alpha ==1) return;
    if (alpha==0){
	for (size_t i=0; i<len; ++i) x[i] = 0;
    } else {
	for (size_t i=0; i<len; ++i) x[i] *= alpha;
    }
}


/** 
 * @brief Skalarprodukt res = x^T * y
 *
 * @param x double Vektor
 * @param y double Vektor
 * @param length Laenge der Vektoren
 * @retrun skalarprodukt der Vektoren
 */
double
blasl1_ddot(const double * x, const double * y, size_t len){
    double res=0;
    for (size_t i=0; i<len; ++i) res += x[i]*y[i];
    return res;
}

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
blasl1_daxpy(double* x, double* y, index length, double alpha, double beta){
    //x = beta*x + alpha*y
    if(beta!=1){
	for(size_t i=0; i<length; i++){
	    x[i] = beta*x[i];
	}
    }

    if(alpha != 0){
	for(size_t i=0; i<length; i++){
	    x[i] += alpha*y[i];
	}
    }
}
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
blasl1_dcopy(const double* a, double* b, index length, double alpha){
    //Kopiert alpha*a in Vektor b
    for(size_t i=0; i<length; i++){
	b[i] = alpha* a[i];
    }
}

/**
 * @brief Index Vektor kopieren b = alpha * a
 * 
 * @param a Vektor
 * @param b Vektor
 * @param length Länge von a und b
 * @param alpha 
 */
void
blasl1_icopy(const index* a, index* b, index length, double alpha){
    //Kopiert alpha*a in Vektor b
    for(size_t i=0; i<length; i++){
	b[i] = alpha* a[i];
    }
}

void blasl1_intcopy(const int* a, int* b, index length, double alpha){
    //Kopiert alpha*a in Vektor b
    for(size_t i=0; i<length; i++){
	b[i] = alpha* a[i];
    }
}

// print functions

/**
 * @brief Print Double-Vektor, jeder Wert in 1 Zeile
 * 
 * @param data Vektor
 * @param len Länge des Vektors
 */
void blasl1_dprint(double * data, size_t len){
    for(size_t i=0;i<len;++i){
	printf("%4.2lf\n",data[i]);
    }
}

/**
 * @brief Print index-Vektor, jeder Wert in 1 Zeile
 * 
 * @param data Vektor
 * @param len Länge des Vektors
 */
void blasl1_iprint(index* data, size_t len){
    for(size_t i=0;i<len;++i){
	printf("%g\n",(double) data[i]);
    }
}

void blasl1_intprint(index* data, size_t len){
    for(size_t i=0;i<len;++i){
	printf("%i\n", data[i]);
    }
}


#endif //BLAS_LEVEL1_H
