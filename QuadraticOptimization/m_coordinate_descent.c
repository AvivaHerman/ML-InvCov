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
#include "m_coordinate_descent.h"


static void M_CD_iteration(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, 
						   double lambda, int index, double* maxDelta, size_t* C, size_t mc);

static double M_CD_selective(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, double lambda, size_t* C, size_t mc);

static double M_debiasing_CD_selective(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, size_t* C, size_t mc);

static void getCoefficient(double* v, double* Mx, double* Mv, double* c, size_t M_size, size_t* C, size_t mc, double* a1, double* a2);

static void update(double* x, double* Mx, double* v, double* Mv, double* c, double alpha, size_t* C, size_t mc, size_t M_size, double lambda);

static void update_Mx(size_t prev_lev_size, size_t curr_lev_size, size_t* C, size_t M_size, double* x, double* Mx, double* M);


double M_CD(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, double lambda)
{
    double maxDelta = 0, a1, a2, alpha, Jopt;
    int i;
	
    for (i = 0; i < M_size; ++i) {
		M_CD_iteration(M, x, Mx, c, v, Mv, M_size, lambda, i, &maxDelta, NULL, M_size);
    }
	
	getCoefficient(v, Mx, Mv, c, M_size, NULL, M_size, &a1, &a2);
	alpha = ApplyLinesearch(0, a1, a2, lambda, 0, 10, x, v, NULL, M_size , &Jopt);
	update(x, Mx, v, Mv, c, alpha, NULL, M_size, M_size, lambda);
	
    return maxDelta;
}


double M_MLCD(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, size_t* C, double lambda, double tol, int coarseningRatio, size_t* level_sizes)
{
    double tmp, maxDelta = 0;
    size_t i, j, wanted_coarse, minCoarseVars = 32, sup_size = 0, part_M_size;
    int lev;
	double *part_M = NULL, *part_x = NULL, *part_c = NULL, *part_Mx = NULL;

	part_Mx = (double *) malloc(sizeof(double) * M_size);
	if (part_Mx == NULL) {
		mexPrintf("malloc failed: part_Mx\n");
		goto out;
	}

    for (i = 0; i < M_size; ++i) {
        if (x[i] != 0) {
            part_Mx[i] = BIG_NUM * (fabs(c[i] - Mx[i]) + fabs(x[i]));
            sup_size++;
        } else {
            part_Mx[i] = fabs(fabs(c[i] - Mx[i]));
        }
    }

    for(i = 0; i < M_size; ++i) {
        C[i] = i;
    }

    quickSort_descend(part_Mx, C, M_size);
	
	free(part_Mx);
	part_Mx = NULL;

    for (i = 0; i < sup_size; ++i) { // beyond this point all xi's should be zero.
        if (x[C[i]] == 0) {
            mexPrintf("WARNING: Index among firsts is not in the support!!!!\n");
        }
    }

    lev = 0;
    wanted_coarse = M_size / coarseningRatio;

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
	
	part_M_size = level_sizes[lev];
	
	part_M = (double *) malloc(sizeof(double) * part_M_size * part_M_size);
	if (part_M == NULL) {
		mexPrintf("malloc failed: part_M\n");
		goto out;
	}
	
	part_x = (double *) malloc(sizeof(double) * part_M_size);
	if (part_x == NULL) {
		mexPrintf("malloc failed: part_x\n");
		goto malloc_part_x_err;
	}
	
	part_Mx = (double *) malloc(sizeof(double) * part_M_size);
	if (part_Mx == NULL) {
		mexPrintf("malloc failed: part_Mx\n");
		goto malloc_part_Mx_err;		
	}
	
	part_c = (double *) malloc(sizeof(double) * part_M_size);
	if (part_c == NULL) {
		mexPrintf("malloc failed: part_c\n");
		goto malloc_part_c_err;		
	}
	
	for (i = 0; i < part_M_size; ++i) {
		for (j = 0; j < part_M_size; ++j) {
			part_M[i * part_M_size + j] = M[C[i] * M_size + C[j]];
		}
		part_x[i] 		= x[C[i]];
		part_Mx[i] 		= Mx[C[i]];
		part_c[i]       = c[C[i]];
	}

	//The exact solution for the partial problem
    maxDelta = M_CD(part_M, part_x, part_Mx, part_c, v, Mv, part_M_size, lambda);
    tmp = M_CD(part_M, part_x, part_Mx, part_c, v, Mv, part_M_size, lambda);
	
    while ((tmp > 0.01 * maxDelta) && (tmp > tol)) {
        tmp = M_CD(part_M, part_x, part_Mx, part_c, v, Mv, part_M_size, lambda);
		
    }

	for (i = 0; i < part_M_size; ++i) {
		x[C[i]]       = part_x[i];
		Mx[C[i]] 	  = part_Mx[i];
	}
	
	free(part_c);
	free(part_Mx);
	free(part_x);
	free(part_M);
    // *******************************************************************
    // *******************************************************************
    
	// The way up:
    for (--lev; lev >= 0; --lev) {
		update_Mx(level_sizes[lev + 1], level_sizes[lev], C, M_size, x, Mx, M);

		quickSort_ascend(NULL, C, level_sizes[lev]);// to make the data more cache friendly
		
        tmp = M_CD_selective(M, x, Mx,c, v, Mv, M_size, lambda, C, level_sizes[lev]);
        maxDelta = tmp;
        while ((tmp > 0.5 * maxDelta) && (tmp > tol)) {
			tmp = M_CD_selective(M, x, Mx, c, v, Mv, M_size, lambda, C, level_sizes[lev]);
        }
    }
    // *******************************************************************
    
	// final level
	update_Mx(level_sizes[0], M_size, C, M_size, x, Mx, M);

    return M_CD(M, x, Mx, c, v, Mv, M_size, lambda);

malloc_part_c_err:
	free(part_Mx);
malloc_part_Mx_err:
	free(part_x);
malloc_part_x_err:
	free(part_M);
out:
	return -1;
}


static void update_Mx(size_t prev_lev_size, size_t curr_lev_size, size_t* C, size_t M_size, double* x, double* Mx, double* M)
{
	size_t i, j;
	double tmp;
	
	for (j = prev_lev_size; j < curr_lev_size; ++j) {
		Mx[C[j]] = 0;
	}
	
	for (i = 0; i < prev_lev_size; ++i) {
		tmp = x[C[i]];
		if (fabs(tmp) > SMALL_NUM) {
			for (j = prev_lev_size; j < curr_lev_size; ++j) {
				Mx[C[j]] += M[C[i] * M_size + C[j]] * tmp; 
			}
		}
	}
}


static double M_CD_selective(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, double lambda, size_t* C, size_t mc)
{
    double maxDelta = 0, a1, a2, alpha, Jopt;
    size_t i;
	
    for(i = 0; i < mc; ++i) {
		M_CD_iteration(M, x, Mx, c, v, Mv, M_size, lambda, C[i], &maxDelta, C, mc);
    }
	
	getCoefficient(v, Mx, Mv, c, M_size, C, mc, &a1, &a2);
	alpha = ApplyLinesearch(0, a1, a2, lambda, 0, 10, x, v, NULL, M_size , &Jopt);
	update(x, Mx, v, Mv, c, alpha, C, mc, M_size, lambda);
	
    return maxDelta;
}

static double M_debiasing_CD_selective(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, size_t* C, size_t mc)
{
    return M_CD_selective(M, x, Mx, c, v, Mv, M_size, 0, C, mc);
}


static void M_CD_iteration(double* M, double* x, double* Mx, double* c, double* v, double* Mv, size_t M_size, 
						   double lambda, int index, double* maxDelta, size_t* C, size_t mc)
{
	double delta, x_old, x_new = 0, Mii, tmp;
	size_t j, offset;
	
	x_old = x[index];
	offset = index * M_size;
	Mii = M[offset + index];

	x_new = shrinkage(lambda / Mii, ((c[index] - Mx[index] - Mv[index]) / Mii) + x_old);
	
	delta = x_new - x_old;
	
	v[index] = delta;

	if (fabs(delta) > SMALL_NUM) {
		*maxDelta = max(fabs(delta), *maxDelta);

		if (C == NULL) {
			for (j = 0; j < M_size; ++j) {
				tmp = (M[offset + j] * delta);
				Mv[j] += tmp;
			}
		} else {
			for (j = 0; j < mc; ++j) {
				tmp = (M[offset + C[j]] * delta);
				Mv[C[j]] += tmp;
			}
		}
	}
	
}


static void getCoefficient(double* v, 			// IN
					       double* Mx, 			// IN
						   double* Mv, 			// IN
					   	   double* c, 			// IN
						   size_t M_size, 		// IN
						   size_t* C, 			// IN
						   size_t mc,			// IN
						   double* a1, 			// OUT
						   double* a2)			// OUT
{
	size_t i;
	*a1 = 0;
	*a2 = 0;
	
	if (C == NULL) {
		for (i = 0; i < M_size; ++i) {
			*a1 += v[i] * (Mx[i] - c[i]);
			*a2 += v[i] * Mv[i];
		}
	} else {
		for (i = 0; i < mc; ++i) {
			*a1 += v[C[i]] * (Mx[C[i]] - c[C[i]]);
			*a2 += v[C[i]] * Mv[C[i]];
		}
	}
	
	*a2 *= 0.5;
}


static void update(double* x, double* Mx, double* v, double* Mv, double* c, double alpha, size_t* C, size_t mc, size_t M_size, double lambda)
{
	size_t i;
	
	if (C == NULL) {
		for (i = 0; i < M_size; ++i) {
			x[i]  	  += alpha * v[i];
			Mx[i]     += alpha * Mv[i];
			v[i]   	   = 0;
			Mv[i]      = 0;
		}
	} else {
		for (i = 0; i < mc; ++i) {
			x[C[i]]  	 += alpha * v[C[i]];
			Mx[C[i]]     += alpha * Mv[C[i]];
			v[C[i]]       = 0;
			Mv[C[i]]  	  = 0;
		}
	}
}
