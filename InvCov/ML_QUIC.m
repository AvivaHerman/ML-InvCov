function [X, W, opt, time, iter, support, numActive, vcycle] = ML_QUIC(mode, S, L, tol, msg, maxIter, X0, W0)
% [X W opt time iter support numActive vcycle] = ML_QUIC("default", S, L, tol, msg, maxIter, X0, W0)
% [X W opt time iter support numActive vcycle] = ML_QUIC("trace", S, L, tol, msg, maxIter, X0, W0)
%
% Arguments:
% mode      selects computation
% S         empirical covariance matrix
% L         matrix of the regularization parameters
% tol       convergence threshold
% msg       verbosity of print statistics
% maxIter   maximum number of Newton iterations to execute
% X0        initial value for the inverse covariance matrix
% W0        initial value for the covariance matrix
%
% Returned values:
% X         inverse covariance matrix, in path mode a vector of matrices
%           the solution of
%             min_X f(X), where f(X) = -logdet(X) + trace(XS) + |Lambda.*X|
% W         inverse of X, in path mode a vector of matrices
% opt       value of f(X) at optimum found; an array in "trace" mode
% time      cpu time; an array in "trace" mode
% iter      the number of Newton iterations executed in "default" or
%           "trace" mode, an array of values in "path" mode
% support   size of the support set; an array in "trace" mode
% numActive size of the active set; an array in "trace" mode
% vcycle    the number of newton iteration in which begins a new vcycle
%
% The "trace" mode allows the computation of approximations to the optimal
% value and the cputime used to acquire them along the way.  Useful for
% plotting graphs showing the speed of convergence.
    
disp('Please compile the executable with make!');

                                         
