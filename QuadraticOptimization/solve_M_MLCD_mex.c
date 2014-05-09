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
#include "mex.h"
#include "m_coordinate_descent.h"
/* mex solve_M_MLCD_mex.c */


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int i, j, k, mc, offset, it, maxit = 0;
    double *M, *x, *Mx, *Mx_new, *c, lambda, *cd, *v = NULL, *Mv = NULL, norml1 = 0;
    size_t *C = NULL, *level_sizes = NULL, M_size;
    double *WU, maxDelta = 0.0, tol, coarseningRatio;
	Trace trace;
	
    M_size 			= mxGetM(prhs[0]);
    
	trace.idx = 1;
	
    M 				= mxGetPr(prhs[0]);
    x 				= mxGetPr(prhs[1]);
    Mx 				= mxGetPr(prhs[2]);
	c 				= mxGetPr(prhs[3]);
    lambda 			= *mxGetPr(prhs[4]);
    maxit  			= (int)*mxGetPr(prhs[5]);
    tol 			= *mxGetPr(prhs[6]);
    coarseningRatio = *mxGetPr(prhs[7]);
    
    // tol - accuracy. The algorithm stops when maximal CD update is below tol.
    // coarseningRatio - the ratio that we reduce the dictionary at each level. generally 0.5. 
      
    /*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(M_size, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(M_size, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(maxit * 2, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(maxit * 2, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
	/*Get a pointer to the data space in our newly allocated memory*/
    WU      = mxGetPr(plhs[0]);
    cd      = mxGetPr(plhs[1]);
	Mx_new 	= mxGetPr(plhs[2]);
	trace.fx_trace 	 = mxGetPr(plhs[3]);
	trace.time_trace = mxGetPr(plhs[4]);
	
	trace.time_trace[0] = clock();

	C = (size_t*) malloc(sizeof(size_t) * M_size);
    if (C == NULL) {
        mexPrintf("malloc failed: C\n");
		goto out;
    }
	
	level_sizes = (size_t*) malloc(sizeof(size_t) * MAX_NUM_OF_LEVELS);
    if (level_sizes == NULL) {
        mexPrintf("malloc failed: level_sizes\n");
		goto malloc_level_sizes_err;
    }
	
	v = (double *) malloc(sizeof(double) * M_size);
	if (v == NULL) {
        mexPrintf("malloc failed: v\n");
		goto malloc_v_err;
	}
	
	Mv = (double *) malloc(sizeof(double) * M_size);
	if (Mv == NULL) {
        mexPrintf("malloc failed: Mv\n");
		goto malloc_Mv_err;
	}

    for (i = 0; i < M_size; i++) {
        cd[i] 	   = x[i];
		Mx_new[i]  = Mx[i];
		v[i] 	   = 0;
		Mv[i]	   = 0;
		norml1 	  += fabs(x[i]);
    }

	trace.fx_trace[0] = 0.5 * innerProduct(x, Mx_new, M_size) - innerProduct(x, c, M_size) + lambda * norml1;
	
    for (it = 0; it < maxit; ++it) {
        maxDelta = M_MLCD(M, cd, Mx_new, c, v, Mv, M_size, C, lambda, tol, (size_t)(1 / coarseningRatio), level_sizes, &trace);
		
		trace.fx_trace[trace.idx] = 0;
		for (i = 0; i < M_size; i++) {
			trace.fx_trace[trace.idx] += 0.5 * cd[i] * Mx_new[i] - cd[i] * c[i] + lambda * fabs(cd[i]);
		}
		
		trace.time_trace[trace.idx] = (clock() - trace.time_trace[0]) / CLOCKS_PER_SEC;
		trace.idx++;
		
        if (maxDelta < tol){
            break;
        }   
    }
	
	trace.time_trace[0] = 0;
	*mxGetPr(plhs[5]) = trace.idx;

    *WU = it + 1;

    // for (; it < maxit ; ++it) {
        // maxDelta = m_Debiasing(M, cd, Mx_new, c ,M_size, lambda);
        // if (maxDelta < 0.2 * tol) {
            // break;
        // }   
    // }

    free(Mv);
malloc_Mv_err:
	free(v);
malloc_v_err:
	free(level_sizes);
malloc_level_sizes_err:
	free(C);
out:
	return;
}
