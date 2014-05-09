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
#include "utility.h"


double shrinkage(double mu, double x)
{
    if (x > mu + SMALL_NUM) {
        return (x - mu);
    }
	
    if (-x > mu + SMALL_NUM) {
        return (x + mu);
    }
	
    return 0; 
}


double innerProduct(double *v,double *u, size_t n)
{
    double ans = 0.0;
    size_t k;
    //#pragma omp parallel for shared(v,u,n) private(k) reduction(+: ans) num_threads(4)
    for (k = 0; k < n; k++) {
        ans += u[k] * v[k];
    }
    return ans;
}


double ApplyLinesearch(double a0, double a1, double a2, double mu, double a, double b, double* x, double* v, size_t* subInd, size_t n, double* J_opt)
{
	double alpha_opt = 0, alpha_l, alpha_r, alpha, J;
	double *alpha_critic;
	size_t *indices_critic, num_critic = 0, up_bound_num_critic = 0;
	double sumx = 0, sumv = 0, rhox = 0, rhov = 0, tmp;
	double sign_t, eps = 1e-15, eps_subsection = 1e-7;
	size_t i, debug = 0,j;
	
	if (!USE_LINESEARCH) {
		return 1;
	}
    
    for (i = 0; i < n; i++) {
        up_bound_num_critic += (v[i] != 0.0);
    }
    
  	alpha_critic = (double *) malloc(sizeof(double) * up_bound_num_critic);
    indices_critic = (size_t *) malloc(sizeof(size_t) * up_bound_num_critic);
	
	if (subInd == NULL) {
		for (i = 0; i < n; i++) {
			tmp = (v[i] != 0.0) ? -x[i] / v[i] : b + 1;
			if (tmp < a + eps || tmp > b - eps) { 
				sign_t = sign(x[i] + (a + eps) * v[i]);
				sumx += sign_t * x[i];
				sumv += sign_t * v[i];
			} else { // alpha is critic, in [a,b]
				sign_t = sign(x[i] + (a+eps) * v[i]);
				rhox += sign_t * x[i];
				rhov += sign_t * v[i];
				alpha_critic[num_critic] = tmp;
          
				indices_critic[num_critic] = i;
				num_critic++;
			}
		}
	} else {
		for (i = 0; i < n; i++) {
			tmp = (v[subInd[i]] != 0.0) ? -x[subInd[i]] / v[subInd[i]] : b + 1;
			if (tmp < a + eps || tmp > b - eps) { 
				sign_t = sign(x[subInd[i]] + (a + eps) * v[subInd[i]]);
				sumx += sign_t * x[subInd[i]];
				sumv += sign_t * v[subInd[i]];
			} else { // alpha is critic, in [a,b]
				sign_t = sign(x[subInd[i]] + (a+eps) * v[subInd[i]]);
				rhox += sign_t * x[subInd[i]];
				rhov += sign_t * v[subInd[i]];
				alpha_critic[num_critic] = tmp;
          
				indices_critic[num_critic] = subInd[i];
				num_critic++;
			}
		}
	}
	
	a0 += mu * sumx;
    a1 += mu * sumv;
	if (num_critic == 0) {
        if (debug) {
            mexPrintf("No Critic point\n");
        }
       
		tmp = -a1 / (2 * a2);
        tmp = max(min(tmp, b), a);
        *J_opt = a2 * tmp * tmp + a1 * tmp + a0;
        alpha_opt = tmp;
        if (debug) {
            mexPrintf("%lf,%lf\n", alpha_opt, *J_opt);
        }
	} else {
		// for (j=0; j<num_critic; ++j) {
			// mexPrintf("index=%d, alpha=%lf\n", indices_critic[j], alpha_critic[j]);
		// }
		// free(alpha_critic);
		// free(indices_critic);
		// return 1;
		// quickSort_ascend(alpha_critic, indices_critic, num_critic);
		quickSort_ascend(alpha_critic, indices_critic, num_critic);
		// mexPrintf("\n");
		// for (j=0; j<num_critic; ++j) {
			// mexPrintf("index=%d, alpha=%lf\n", indices_critic[j], alpha_critic[j]);
		// }
		// free(alpha_critic);
		// free(indices_critic);
		// return 1;
		
		
		// for (j=0; j<num_critic; ++j) {
			// mexPrintf("index=%d, alpha=%lf\n", indices_critic[j], alpha_critic[j]);
		// }
        if (debug) {
            mexPrintf("num_critic = %d\n", num_critic);
            mexPrintf("a = %lf < ", a);
            for (i = 0; i < num_critic; i++) {
                mexPrintf("%2.14lf < ", alpha_critic[i]);
            }
            mexPrintf("%2.14lf\n", b);
        }
		alpha_l = a;
		alpha_r = alpha_critic[0];
        if (debug) {
            mexPrintf("Checking out section (%2.14lf,%2.14lf)\n", alpha_l, alpha_r);
        }
		tmp = (mu * rhov + a1);
		alpha_opt = - tmp / (2 * a2);
		alpha_opt = max(min(alpha_r, alpha_opt), alpha_l);
		*J_opt = (a2 * alpha_opt + tmp) * alpha_opt + a0 + mu * rhox;
        if (debug) {
            mexPrintf("(alpha_opt,J_opt) = (%2.14lf,%2.14lf)\n", alpha_opt, *J_opt);
        }
		sign_t = sign(x[indices_critic[0]] + (a + eps) * v[indices_critic[0]]);
		rhox = rhox - 2 * x[indices_critic[0]] * sign_t;
        rhov = rhov - 2 * v[indices_critic[0]] * sign_t;
		i = 1;
		
		while (fabs(alpha_critic[i] - alpha_critic[i - 1]) < eps_subsection && i < num_critic) {
            sign_t = sign(x[indices_critic[i]] + (a + eps) * v[indices_critic[i]]);
			rhox = rhox - 2 * x[indices_critic[i]] * sign_t;
			rhov = rhov - 2 * v[indices_critic[i]] * sign_t;
            i++;
		}
		
		while (i < num_critic) {
			alpha_l = alpha_critic[i - 1];
			alpha_r = alpha_critic[i];
            if (debug) {
                mexPrintf("Checking out section (%2.14lf,%2.14lf)\n", alpha_l, alpha_r);
            }
			tmp = (mu * rhov + a1);
			alpha = - tmp / (2 * a2);
			alpha = max(min(alpha_r, alpha), alpha_l);
			J = (a2 * alpha + tmp) * alpha + a0 + mu * rhox;
            if (debug) {
                mexPrintf("(alpha,J) = (%2.14lf,%2.14lf)\n", alpha, J);
            }
			if ( J <= *J_opt) {
				*J_opt = J;
				alpha_opt = alpha;
			} else { // optional? else break. A shortcut because of the convexity.
                break;
            }
			sign_t = sign(x[indices_critic[i]] + (a + eps) * v[indices_critic[i]]);
			rhox = rhox - 2 * x[indices_critic[i]] * sign_t;
			rhov = rhov - 2 * v[indices_critic[i]] * sign_t;
            i++;
			while (fabs(alpha_critic[i] - alpha_critic[i - 1]) < eps_subsection && i < num_critic) {
				sign_t = sign(x[indices_critic[i]] + (a + eps) * v[indices_critic[i]]);
                rhox = rhox - 2 * x[indices_critic[i]] * sign_t;
				rhov = rhov - 2 * v[indices_critic[i]] * sign_t;
                i++;
			}
		}
		
		alpha_l = alpha_critic[num_critic - 1];
		alpha_r = b;
		tmp = (mu * rhov + a1);
		alpha = - tmp / (2 * a2);
		alpha = max(min(alpha_r, alpha), alpha_l);
		J = (a2 * alpha + tmp) * alpha + a0 + mu * rhox;
		if ( J < *J_opt) {
			*J_opt = J;
			alpha_opt = alpha;
		}
	}
	free(alpha_critic);
	free(indices_critic);
	return alpha_opt;
}


double sign(double a)
{
    return 2 * (a >= 0) - 1; // this returns 1 for 0 as well...
}
