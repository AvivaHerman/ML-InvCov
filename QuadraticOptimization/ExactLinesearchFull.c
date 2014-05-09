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
#include "utility.h"
/* ExactLinesearchFull.c */


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // input parameters:
    double a2, a1, a0, a, b, mu;
    double *x, *v;
    int n_xv;
	
    // output parameters:
    double *alpha_opt, *J_opt;
    // help variables
   
    x = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    
    a2 = *mxGetPr(prhs[2]);
    a1 = *mxGetPr(prhs[3]);
    a0 = *mxGetPr(prhs[4]);
    a  = *mxGetPr(prhs[5]);
    b  = *mxGetPr(prhs[6]);
    mu = *mxGetPr(prhs[7]);
	
	//Allocate memory and assign output pointer
    plhs[0]   = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1]   = mxCreateDoubleMatrix(1, 1, mxREAL);
    alpha_opt = mxGetPr(plhs[0]);
    J_opt     = mxGetPr(plhs[1]);
    
    if (mxGetM(prhs[0]) != mxGetM(prhs[1])) {
        printf("error: x and v are not from the same size! Aborting!\n");
        return;
    }
	
    if (mxGetN(prhs[0]) != mxGetN(prhs[1])) {
        printf("error: x and v are not from the same size! Aborting!\n");
        return;
    }
	
    n_xv = mxGetM(prhs[0]) * mxGetN(prhs[0]);
    // END OF MEX INTERFACE ISSUES.
    *alpha_opt = ApplyLinesearch(a0, a1, a2,  mu,  a,  b,  x,  v, 0, n_xv, J_opt);
}
