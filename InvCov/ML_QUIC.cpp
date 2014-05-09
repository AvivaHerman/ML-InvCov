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

// This is the MEX wrapper for ML_QUIC.  The algorithm is in QUIC_utils.C.

// Invocation form within Matlab or Octave:
// [X W opt time iter] = QUIC(mode, ...)
// [X W opt time iter] = QUIC("default", S, L, tol, msg, maxIter,
//                            X0, W0)
// [X W opt time iter] = QUIC("path", S, L, path, tol, msg, maxIter,
//                            X0, W0)
// [X W opt time iter] = QUIC("trace", S, L, tol, msg, maxIter,
//                                 X0, W0)
// See the README file and QUIC.m for more information.

#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "QUIC_utils.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("QUIC:arguments",
                "Missing arguments, please specify\n"
                "             S - the empirical covariance matrix, and\n"
                "             L - the regularization parameter.");
    }

	fX_Info *fx_info_ptr;
	fX_Info fx_info;
	fX_InfoFullLambda fx_info_lambda;
	
	fx_info_ptr = &fx_info_lambda;
	
	Trace_Info trace;
	
    long argIdx = 0;
    
	char mode[8];
    mode[0] = 'd';
    if (mxIsChar(prhs[0])) {
        mxGetString(prhs[0], mode, 8);
        if (strcmp(mode, "path") &&
                strcmp(mode, "trace") &&
                strcmp(mode, "default"))
            mexErrMsgIdAndTxt("QUIC:arguments",
                    "Invalid mode, use: 'default', 'path' or "
                    "'trace'.");
        argIdx++;
    }
    
	// The empirical covariance matrix:
    if (!mxIsDouble(prhs[argIdx]))
        mexErrMsgIdAndTxt("QUIC:type",
                "Expected a double matrix. (Arg. %d)",
                argIdx + 1);
    const double* S = mxGetPr(prhs[argIdx]);
    uint32_t p = (uint32_t) mxGetN(prhs[argIdx]);
    if (p != uint32_t(mxGetM(prhs[argIdx]))) {
        mexErrMsgIdAndTxt("QUIC:dimensions",
                "Expected a square empirical covariance matrix.");
    }
    argIdx++;
    
    // Regularization parameter matrix:
    if (!mxIsDouble(prhs[argIdx]))
        mexErrMsgIdAndTxt("QUIC:type",
                "Expected a double matrix. (Arg. %d)",
                argIdx + 1);

	if (mxGetN(prhs[argIdx]) == 1 && mxGetM(prhs[argIdx]) == 1) {
		fx_info_ptr = &fx_info;
		fx_info_ptr->Lambda = mxGetPr(prhs[argIdx]);	
    } else {
        if (mxGetN(prhs[argIdx]) != p && mxGetM(prhs[argIdx]) != p) {
            mexErrMsgIdAndTxt("QUIC:dimensions",
                    "The regularization parameter is not a scalar\n"
                    "              or a matching matrix.");
        }
		fx_info_ptr = &fx_info_lambda;
        fx_info_ptr->Lambda = mxGetPr(prhs[argIdx]);
    }
    argIdx++;
    
    double tol = 1e-6;
    if (nrhs > argIdx) {
        if (!mxIsDouble(prhs[argIdx]))
            mexErrMsgIdAndTxt("QUIC:type",
                    "Expected a double scalar. (Arg. %d)",
                    argIdx + 1);
        tol = mxGetScalar(prhs[argIdx]);
        if (tol < 0) {
            mexErrMsgIdAndTxt("QUIC:tol",
                    "Negative tolerance value.");
        }
        argIdx++;
    }
    
    uint32_t msg = QUIC_MSG_FAILURE;
    if (nrhs > argIdx) {
        msg = (uint32_t) mxGetScalar(prhs[argIdx]);
        argIdx++;
    }
    
    // Maximum number of Newton ierations (whole matrix update):
    uint32_t maxIter = 1000;
    if (nrhs > argIdx) {
        maxIter = (uint32_t) mxGetScalar(prhs[argIdx]);
        argIdx++;
    }
    
    double* X0 = NULL;
    double* W0 = NULL;
    if (nrhs > argIdx) {
        if (!mxIsDouble(prhs[argIdx]))
            mexErrMsgIdAndTxt("QUIC:type",
                    "Expected a double matrix. (Arg. %d)",
                    argIdx + 1);
        if (p != mxGetM(prhs[argIdx]) || p != mxGetN(prhs[argIdx]))
            mexErrMsgIdAndTxt("QUIC:dimensions",
                    "Matrix dimensions should match.\n"
                    "             Maybe incorrect mode is specified?");
        X0 = mxGetPr(prhs[argIdx]);
        argIdx++;
        if (nrhs == argIdx)
            mexErrMsgIdAndTxt("QUIC:initializations",
                    "Please specify both the initial estimate\n"
                    "             and the inverse.\n"
                    "             Maybe incorrect mode is specified?");
        if (!mxIsDouble(prhs[argIdx]))
            mexErrMsgIdAndTxt("QUIC:type",
                    "Expected a double matrix. (Arg. %d)",
                    argIdx + 1);
        if (p != mxGetM(prhs[argIdx]) || p != mxGetN(prhs[argIdx])) {
            mexErrMsgIdAndTxt("QUIC:dimensions",
                    "Matrix dimensions should match.\n"
                    "             Maybe incorrect mode is specified?");
        }
        W0 = mxGetPr(prhs[argIdx]);
        argIdx++;
    }
    
    double* X = NULL;
    double* W = NULL;

	mwSize dims[] = {p, p};
	mxArray* tmp = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
	X = (double *) mxGetPr(tmp);	
	if (!X) {
		MSG("Failure: X is NULL\n");
		return;
	}
	if (nlhs > 0)
		plhs[0] = tmp;
		
	tmp = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
	W = (double *) mxGetPr(tmp);
	if (!W) {
		MSG("Failure: W is NULL\n");
		return;
	}
	if (nlhs > 1)
		plhs[1] = tmp;

    if (X0 != NULL) {
        memcpy(X, X0, sizeof(double) * p * p);
        memcpy(W, W0, sizeof(double) * p * p);
    } else {
        memset(X, 0, sizeof(double) * p * p);
        memset(W, 0, sizeof(double) * p * p);
        for (unsigned long i = 0; i < p * p; i += (p + 1)) {
            X[i] = 1.0;
            W[i] = 1.0;
        }
    }
	
    double* opt				= NULL;
    double* cputime 		= NULL;
    uint32_t* iter 			= NULL;
	uint32_t* suppArr 		= NULL;
	uint32_t* activeNumArr 	= NULL;
	uint32_t* vcycleIter 	= NULL;
	
    unsigned long traceSize = VSYCLE_LEVELS;
    unsigned long iterSize = 1;
	
    if (mode[0] == 't')
        traceSize = maxIter * VSYCLE_LEVELS;
		
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(traceSize, 1, mxREAL);
        opt = mxGetPr(plhs[2]);
    }
	
    if (nlhs > 3) {
        plhs[3] = mxCreateDoubleMatrix(traceSize, 1, mxREAL);
        cputime = mxGetPr(plhs[3]);
    }
	
    if (nlhs > 4) {
        mwSize dims[] = {iterSize};
        plhs[4] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
        iter = (uint32_t *) mxGetData(plhs[4]);
    }
	
	if (nlhs > 5) {
        mwSize dims[] = {traceSize};
        plhs[5] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
        suppArr = (uint32_t *) mxGetData(plhs[5]);
    }
	
	if (nlhs > 6) {
        mwSize dims[] = {traceSize};
        plhs[6] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
        activeNumArr = (uint32_t *) mxGetData(plhs[6]);
    }
	
	if (nlhs > 7) {
        mwSize dims[] = {traceSize};
        plhs[7] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
        vcycleIter = (uint32_t *) mxGetData(plhs[7]);
    }
	
	fx_info_ptr->p 		 = p;
	fx_info_ptr->l1normX = 0.0;
	fx_info_ptr->trSX 	 = 0.0;
	fx_info_ptr->logdetX = 0.0;
	fx_info_ptr->fX 	 = INITIAL_FX;
	fx_info_ptr->X 		 = X;
	fx_info_ptr->S 		 = S;
	fx_info_ptr->W		 = W;
	fx_info_ptr->tol	 = tol;
	
	trace.opt 		   = opt;
	trace.iter 		   = iter;
	trace.vcycleIter   = vcycleIter;
	trace.cputime 	   = cputime;
	trace.suppArr 	   = suppArr;
	trace.activeNumArr = activeNumArr;
	trace.traceIdx 	   = 0;
	trace.mode 		   = mode[0];
	trace.msg 		   = msg;

    QUIC(*fx_info_ptr, trace, maxIter, ML_QUIC);
}
