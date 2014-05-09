function [A,X] = generateRandomPlanarGraphJumpy(N)
[A,X] = generateRandomPlanarGraph(N);
C = generateDeCompositionMatrix(A);
n_edges = size(C,1);
W = rand(n_edges,1);
W = spdiags(W,0,n_edges,n_edges);
A = C'*W*C;
return;

function C = generateDeCompositionMatrix(A)
[I,J,~] = find(A);
S = (I>J)';
I = I(S)';
J = J(S)';
D = [ones(size(I)),-ones(size(I))];
C = sparse([1:length(I),1:length(I)],[I,J],D);
return;