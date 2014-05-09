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
#include "coordinate_descent.h"
/* mex solveL2L1_MLCD_mex.c */


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int i, j, k, mc, offset, it, maxit = 0;
    double *A, *x, *r, *r_new, mu, *cd, *Atr_prev, *Av = NULL, *v = NULL, *prev_r = NULL, norml1 = 0;
    size_t *C = NULL, *level_sizes = NULL;
    double *WU, maxDelta = 0, tol, *Atr, coarseningRatio;
	Trace trace;
	
    size_t rowLen = mxGetN(prhs[0]);
    size_t colLen = mxGetM(prhs[0]);
    
	trace.idx = 1;
	
    A 				= mxGetPr(prhs[0]);
    x 				= mxGetPr(prhs[1]);
    r 				= mxGetPr(prhs[2]);
    mu 				= *mxGetPr(prhs[3]);
    maxit 			= (int)*mxGetPr(prhs[4]);
    tol   			= *mxGetPr(prhs[5]);
    coarseningRatio = *mxGetPr(prhs[6]);
    Atr_prev 		= mxGetPr(prhs[7]);
    
    // mu is sometimes called lambda
    // tol - accuracy. The algorithm stops when maximal CD update is below tol.
    // coarseningRatio - the ration that we reduce the dictionary at each level. generally 0.5. 
	
    /*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rowLen, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(rowLen, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(maxit * 2, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(maxit * 2, 1, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/*Get a pointer to the data space in our newly allocated memory*/
    WU    			 = mxGetPr(plhs[0]);
    cd    			 = mxGetPr(plhs[1]);
    Atr    			 = mxGetPr(plhs[2]);
    r_new   		 = mxGetPr(plhs[3]);
	trace.fx_trace 	 = mxGetPr(plhs[4]);
	trace.time_trace = mxGetPr(plhs[5]);
	
	trace.time_trace[0] = clock();
	
	C = (size_t *) malloc(sizeof(size_t) * rowLen);
    if (C == NULL) {
        mexPrintf("malloc failed: C\n");
		goto out;
    }
    
	level_sizes = (size_t *) malloc(sizeof(size_t) * MAX_NUM_OF_LEVELS);
	if (level_sizes == NULL) {
        mexPrintf("malloc failed: level_sizes\n");
		goto malloc_level_sizes_err;
    }
	
	v = (double *) malloc(sizeof(double) * rowLen);
	if (v == NULL) {
        mexPrintf("malloc failed: v\n");
		goto malloc_v_err;
    }
	
	Av = (double *) malloc(sizeof(double) * colLen);
	if (Av == NULL) {
        mexPrintf("malloc failed: Av\n");
		goto malloc_Av_err;
    }
	
	prev_r = (double *) malloc(sizeof(double) * colLen);
	if (prev_r == NULL) {
        mexPrintf("malloc failed: prev_r\n");
		goto malloc_prev_r_err;
    }
    
    for (i = 0; i < colLen; i++) {
        r_new[i]  = r[i];
		prev_r[i] = r[i];
		Atr[i]    = Atr_prev[i];
		Av[i]     = 0;
    }
	
    for (i = 0; i < rowLen; i++) {
        cd[i]   = x[i];
		v[i]    = 0;
		norml1 += fabs(x[i]);
    }
	
	trace.fx_trace[0] = 0.5 * innerProduct(r, r, colLen) + mu * norml1;
	
	mexPrintf("MLCD: %lf\n", trace.fx_trace[0]);
	
    for (it = 0; it < maxit; ++it) {
        maxDelta = MLCD(A, cd, v, Av, r_new, prev_r, Atr, C, rowLen, colLen, mu, tol, (size_t)(1 / coarseningRatio), level_sizes, &trace);
		
		norml1 = 0;
		for (i = 0; i < rowLen; i++) {
			norml1 += fabs(cd[i]);
		}
		
		trace.fx_trace[trace.idx] = 0.5 * innerProduct(r_new, r_new, colLen) + mu * norml1;
		trace.time_trace[trace.idx] = (clock() - trace.time_trace[0]) / CLOCKS_PER_SEC;
		trace.idx++;
		
        if (maxDelta < tol) {
            break;
        }   
    }
	
	trace.time_trace[0] = 0;
	*mxGetPr(plhs[6]) = trace.idx;

    *WU = it + 1;
	
    // for (; it < maxit ; ++it) {
        // maxDelta = Debiasing(A, cd, Av, Atr, rowLen, colLen, mu);
        // if (maxDelta < 0.2*tol) {
            // break;
        // }   
    // }
	
	free(prev_r);
malloc_prev_r_err:
	free(Av);
malloc_Av_err:
	free(v);
malloc_v_err:
    free(level_sizes);
malloc_level_sizes_err:
	free(C);
out:
	return;
}
