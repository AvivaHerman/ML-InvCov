invocation of multilevel coordinate descent algorithm for L1 panalized least-squares:

This the README file for the Matlab release of the solve_M_CD_mex, solve_M_MLCD_mex programs.

To compile the Matlab executable:

> make

which produces a 64-bit MEX file.

The arguments and return values when invoked from within Matlab are
documented in solve_M_MLCD_mex.m and solve_M_CD_mex.m. 
You may type 'help solve_M_MLCD_mex' or 'help solve_M_CD_mex' to obtain the full
documentation.

Sample use:
> [x, Mx_res, iter] = solve_M_CD_mex(M, x, Mx, c, lambda, max_iter, accuracy);
> [iter, x_res, Mx_res] = solve_M_MLCD_mex(M, x, Mx, c, lambda, max_iter, accuracy, ratio);