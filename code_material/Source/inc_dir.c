#include "hpc.h"

void inc_dir_r(double* r, const index* dir_ind, const size_t n_dir){
	for (size_t i = 0; i < n_dir; ++i){
		r[dir_ind[i]] = 0;
	}
}

void inc_dir_u(double* u, const double* dir, const index* dir_ind, const size_t n_dir){
	for (size_t i = 0; i < n_dir; ++i){
		u[dir_ind[i]] = dir[i];
	}
}