function [ A_true,X ] = generateExperiment( n,k )
addpath('GraphsLap');

% A_true = speye(n);
% A_true = generateRandomPlanarGraph(n);

A_true = delsq(numgrid('C',ceil(sqrt(n)))) ;
n = size(A_true,2);

% A_true = sprandn(n,n,4/n,1/n);
% A_true = A_true'*A_true;
% n = size(A_true,2);


A_true = A_true + 1e-8*speye(size(A_true));
% rescaling = spdiags(1./sqrt(diag(A_true)),0,n,n);
% A_true = (rescaling*A_true*rescaling);

% generating X:
L_true = chol(A_true);
Y = randn(n,k);
mu_hat = sum(Y,2)/k;
Y = Y - mu_hat*ones(1,k);

X = L_true\Y; % we need to multiply Y by the sqrt(covariance) and A_true is the inverse of the cov.
% norm((1/k)*X*X' - A_true,'fro')/(n) % this is not true for a partial rank
% X

% norm(((1/k)*X*X' - A_true).*(A_true~=0),'fro')/(n) % this is not true for a partial rank
% X
rmpath('GraphsLap');
end

