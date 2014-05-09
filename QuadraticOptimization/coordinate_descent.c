/*---------------------------------------------------------------------------
Copyright (2014): Eran Treister, Aviva Herman and Irad Yavneh. 
This code is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."
---------------------------------------------------------------------------*/
#include "coordinate_descent.h"


static void CD_iteration(double* A, double* x, double* v, double* Av, double* r_new, double* Atr, size_t colLen, double mu, int index, double* maxDelta);

static double CD_selective(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t* C, size_t mc, size_t colLen, size_t rowlen, double mu);

static double debiasing_CD_selective(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t* C, size_t mc,  size_t colLen);

static void getCoefficient(double* v, double* prev_r, double* Av, size_t collen, double* a1, double* a2);

static void update(double* A, double* x, double* prev_r, double* r, double* v, double* Av, double alpha, size_t collen, size_t rowlen, double mu);


double CD(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t rowLen,  size_t colLen, double mu)
{
    double maxDelta = 0, a1, a2, alpha,Jopt;
    size_t i;
    
	for (i = 0; i < rowLen; ++i) {
        CD_iteration(A, x, v, Av, r_new, Atr, colLen, mu, i, &maxDelta);
    }
	
	getCoefficient(v, prev_r, Av, colLen, &a1, &a2);
	alpha = ApplyLinesearch(0, a1, a2, mu, 0, 10, x, v, NULL, rowLen , &Jopt);
	update(A, x, prev_r, r_new, v, Av, alpha, colLen, rowLen, mu);
	
    return maxDelta;
}


double MLCD(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t* C, size_t rowLen, size_t colLen, double mu, double tol, int coarseningRatio, size_t* level_sizes)
{
    // THIS IS VCYCLE ONLY!!! 
	// Coarse-level - low level.
    double tmp, maxDelta = 0;
    size_t i, sup_size, wanted_coarse, minCoarseVars = 8;
    int lev;
    sup_size = 0;
	
    for (i = 0; i < rowLen; ++i) {
        if (x[i] != 0) {
            Atr[i] = BIG_NUM * (fabs(Atr[i]) + fabs(x[i]));
            sup_size++;
        } else {
            Atr[i] = fabs(Atr[i]);
        }
    }
	
    for(i = 0; i < rowLen; ++i) {
        C[i] = i;
    }
	
    quickSort_descend(Atr, C, rowLen);
	
    for (i = 0; i < sup_size; ++i) { // beyond this point all xi's should be zero.
        if (x[C[i]] == 0) {
            mexPrintf("WARNING: Index among firsts is not in the support!!!!\n");
        }
    }
	
    lev = 0;
    wanted_coarse = rowLen / coarseningRatio;
	
    // Here we calculate the sizes of the sets: 
    do {
        if (wanted_coarse <= sup_size) {
            level_sizes[lev++] = sup_size;
            break;
        } else {
            level_sizes[lev++] = wanted_coarse;
        }
        wanted_coarse = max(wanted_coarse / coarseningRatio, sup_size);
    } while ((wanted_coarse > minCoarseVars) && (sup_size <= wanted_coarse));
    // The integer values in C from 0 to level_sizes[lev] are the indices of level lev.
    
    quickSort_ascend(NULL, C, level_sizes[--lev]);
	
    // *******************************************************************
    // Exact solve on lowest level:
    // We apply iterations until problem is solved "accurately enough".
    // Stopping condition:
    // 1) MaxDelta on a given iteration < tol.
    // 2) Let MD be the MaxDelta on first iteration. We apply iterations until MaxDelta is below 0.5*MD. 
    maxDelta = CD_selective(A, x, v, Av, r_new, prev_r, Atr, C, level_sizes[lev], colLen, rowLen, mu);
    tmp = CD_selective(A, x, v, Av, r_new, prev_r, Atr, C, level_sizes[lev], colLen, rowLen, mu);
    while ((tmp > 0.5 * maxDelta) && (tmp > tol)) {
        tmp = CD_selective(A, x, v, Av, r_new, prev_r, Atr, C, level_sizes[lev], colLen, rowLen, mu);
    }
    // *******************************************************************
    // *******************************************************************
    
	// The way up:
    for (--lev; lev >= 0; --lev) {
        quickSort_ascend(NULL, C, level_sizes[lev]); // to make the data more cache friendly
        tmp = CD_selective(A, x, v, Av, r_new, prev_r, Atr, C, level_sizes[lev], colLen, rowLen, mu);
        maxDelta = tmp;
        while ((tmp > 0.5 * maxDelta) && (tmp > tol)) {
			tmp = CD_selective(A, x, v, Av, r_new, prev_r, Atr, C, level_sizes[lev], colLen, rowLen, mu);
        }
    }
    // *******************************************************************
    
	// final level
    return CD(A, x, v, Av, r_new, prev_r, Atr, rowLen, colLen, mu);
}


static double CD_selective(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t* C, size_t mc, size_t colLen, size_t rowLen, double mu)
{
    double maxDelta = 0, a1, a2, alpha, Jopt;
    size_t i;
	
    for(i = 0; i < mc; ++i) {
		CD_iteration(A, x, v, Av, r_new, Atr, colLen, mu, C[i], &maxDelta);
    }
	
	getCoefficient(v, prev_r, Av, colLen, &a1, &a2);
	alpha = ApplyLinesearch(0, a1, a2, mu, 0, 10, x, v, NULL, rowLen , &Jopt);
	update(A, x, prev_r, r_new, v, Av, alpha, colLen, rowLen, mu);
	
    return maxDelta;
}


static double debiasing_CD_selective(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t* C, size_t mc, size_t colLen, size_t rowLen)
{
	return CD_selective(A, x, v, Av, r_new, prev_r, Atr, C, mc, colLen, rowLen, 0);
}


static void CD_iteration(double* A, double* x, double* v, double* Av, double* r_new, double* Atr, size_t colLen, double mu, int index, double* maxDelta)
{
	double tmp, x_new = 0, x_old = 0, delta = 0;
	size_t j, offset;
	
	x_old = x[index];
	offset = index * colLen;
	Atr[index] = innerProduct(A + offset, r_new, colLen);
	x_new = shrinkage(mu, Atr[index] + x_old); // here we assume that the dictionary is normalized  A_i^T*A_i = 1

	delta = x_new - x_old;
	
	v[index] = delta;

	if (fabs(delta) > SMALL_NUM) {
		*maxDelta = max(fabs(delta), *maxDelta);
		
		for (j = 0; j < colLen; ++j) {
			tmp 	  = A[offset + j] * delta;
			r_new[j] -= tmp;
			Av[j] 	 -= tmp;
		}
	}
}


static void getCoefficient(double* v, 			// IN
					       double* prev_r, 		// IN
						   double* Av, 			// IN
						   size_t colLen, 		// IN
						   double* a1, 			// OUT
						   double* a2)			// OUT
{
	size_t i;
	
	*a1 = innerProduct(Av, prev_r, colLen);
	*a2 = 0.5 * innerProduct(Av, Av, colLen);
}


static void update(double* A, double* x, double* prev_r, double* r, double* v, double* Av, double alpha, size_t colLen, size_t rowLen, double mu)
{
	size_t i;

	for (i = 0; i < colLen; ++i) {
		r[i]      = prev_r[i] + alpha * Av[i];
		prev_r[i] = r[i];
		Av[i]     = 0;
	}
	
	for (i = 0; i < rowLen; ++i) {
		x[i]   += alpha * v[i];
		v[i]    = 0;
	}
}
