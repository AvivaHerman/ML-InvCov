function [iter, x_res, Mx_res] = solve_M_MLCD_mex(M, x, Mx, c, lambda, max_iter, accuracy, ratio);
%	Let A be the dictionary and y the signal, then:
%	INPUT:
%	M 		 = A' * A (size of m * m)
%	c 		 = A' * y
%	x, Mx 	 = zeroed vectors of size (m * 1)
%	lambda 	 = sparsity regularization parameter
%	max_iter = maximum iterations of Vcycle
%	accuracy = the wanted accuracy
%	ratio	 = the ratio that we reduce the dictionary at each level
%
%	OUTPUT:
%	iter 	 = the number of iteration that were executed until convergence
%	x_res	 = sparse x the minimize F(x)
%	Mx_res	 = M * x_res
   disp('You must run make first!');