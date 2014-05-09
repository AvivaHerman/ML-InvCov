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
#ifndef UTILITY_H
#define UTILITY_H

#include <stddef.h>	// size_t
#include <stdio.h>	// printf
#include <stdlib.h>	// malloc, free, max
#include <math.h>   // fabs
#include <mex.h>	// mexPrintf
#include <time.h>	// clock(), CLOCKS_PER_SEC
#include "sort.h"

#define BIG_NUM 			100000
#define SMALL_NUM 			1e-15
#define MAX_NUM_OF_LEVELS 	10
#define USE_LINESEARCH 		1

typedef struct {
	double* fx_trace;
	double* time_trace;
	int 	idx;
} Trace;

double shrinkage(double mu, double x);

double innerProduct(double *v, double *u, size_t n);

double ApplyLinesearch(double a0, double a1, double a2, double mu, double a, double b, double* x, double* v, int* subInd, int n,double* J_opt);

double sign(double a);

#endif /* UTILITY_H */
