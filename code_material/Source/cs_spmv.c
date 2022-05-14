#include "hpc.h"
#include <stdbool.h>
//y=beta*y + alpha*A*x
index cs_spmv(const cs *A, const double *x, double *y, double alpha,
	double beta){
    index p, j, m, n, nz, *Ap, *Ai ;
    double *Ax, tmp ;

    if (!A || !x || !y) return (0);
    m = A->m; nz = A->nz ; Ap = A->p ; Ai = A->ind ; Ax = A->x ;

    if(beta!=1){
	for(size_t i=0;i<m;++i){
	    y[i]*=beta;
	}
    }

    for (j = 0 ; j < nz ; j++){
	y[Ai[j]]+=alpha*Ax[j]*x[Ap[j]];
    }
    return (1) ;
}

//y=beta*y + alpha*A^T*x
index cs_spmv_transp(const cs *A, const double *x, double *y, double alpha,
	double beta){
    index p, j, m, n, nz, *Ap, *Ai ;
    double *Ax, tmp ;

    if (!A || !x || !y) return (0);
    m = A->m; nz = A->nz ; Ap = A->p ; Ai = A->ind ; Ax = A->x ;

    if(beta!=1){
	for(size_t i=0;i<m;++i){
	    y[i]*=beta;
	}
    }

    for (j = 0 ; j < nz ; j++){
        y[Ap[j]]+=alpha*Ax[j]*x[Ai[j]];
    }
    return (1) ;
}


//Berechnet fuer cs=triplet Format y=beta*y + alpha*A*x
index cs_spmv_const(const cs *A, const double *x, double *y,double alpha,
	double beta, index* fixed, index nfixed)
{
    index p, j, m, n, nz, *Ap, *Ai;
    double *Ax, tmp;
    bool bfixed=false;

    if (!A || !x || !y) return (0);

    m = A->m; nz = A->nz; Ap = A->p; Ai = A->ind; Ax = A->x;

    if(beta!=1){
	for(size_t i=0;i<m;++i){
	    y[i]*=beta;
	}
    }

    if(nfixed<=0){
	return cs_spmv(A,x,y,alpha,beta);
    }else{
	for (j=0; j<nz; j++){
	    for(size_t k=0;k<nfixed;++k){
		if(fixed[k]==Ai[j]){
		    bfixed=true;
		    break;
		}
	    }
	    if(!bfixed){
		y[Ai[j]]+=alpha*Ax[j]*x[Ap[j]];
	    }
	    bfixed=false;
	}
    }
    return (1) ;
}

//Berechnet fuer cs=triplet Format y=beta*y + alpha*A*x fuer A^(T)
index cs_spmv_const_transp(const cs *A, const double *x, double *y,double alpha,
	double beta, index* fixed, index nfixed)
{
    index p, j, m, n, nz, *Ap, *Ai;
    double *Ax, tmp;
    bool bfixed=false;

    if (!A || !x || !y) return (0);

    m = A->m; nz = A->nz; Ap = A->p; Ai = A->ind; Ax = A->x;

    if(beta!=1){
	for(size_t i=0;i<m;++i){
	    y[i]*=beta;
	}
    }

    if(nfixed<=0){
	for (j=0; j<nz; j++){
	    y[Ap[j]]+=alpha*Ax[j]*x[Ai[j]];
	}
    }else{
	for (j=0; j<nz; j++){
	    for(size_t k=0;k<nfixed;++k){
		if(fixed[k]==Ap[j]){
		    bfixed=true;
		    break;
		}
	    }
	    if(!bfixed){
		y[Ap[j]]+=alpha*Ax[j]*x[Ai[j]];
	    }
	    bfixed=false;
	}
    }
    return (1) ;
}
