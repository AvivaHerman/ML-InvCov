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
#ifndef CD_H
#define CD_H


#include "utility.h"
#include "sort.h"

double CD(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t rowLen,  size_t colLen, double mu);

double MLCD(double* A, double* x, double* v, double* Av, double* r_new, double* prev_r, double* Atr, size_t* C, size_t rowLen, size_t colLen, double mu, double tol, int coarseningRatio, size_t* level_sizes);

#endif /* CD_H */
