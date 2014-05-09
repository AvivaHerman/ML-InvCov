function [A,X] = generateRandomCubicGraph(N)
x = rand(N,1);
y = rand(N,1);
z = rand(N,1);
[y,I] = sort(y);
x = x(I);
z = z(I);
X{1} = x;
X{2} = y;
X{3} = z;

TRI = delaunay(x,y,z);

connectionsPerElement = 6;
D = ones(connectionsPerElement*size(TRI,1),3);
for k = 1:size(TRI,1)
    indices = (connectionsPerElement*(k-1)+1) : (connectionsPerElement*k);
    D(indices,1:2) = [TRI(k,1),TRI(k,2) ; TRI(k,2),TRI(k,3) ; TRI(k,1),TRI(k,3) ;...
                      TRI(k,1),TRI(k,4) ; TRI(k,2),TRI(k,4) ; TRI(k,3),TRI(k,4) ];
end

A = sparse(D(:,1),D(:,2),D(:,3),N,N);clear D;

A = A + A';
A = double(A>0);
D = spdiags(sum(A)',0,N,N);
A = -A + D;
% hb = 1/N^(1/3);
% boundary = x < hb | x > 1-hb | y < hb | y > 1-hb | z < hb | z > 1-hb;
% sum(boundary) / N
% A = A(~boundary,~boundary);
% N = size(A,1);
% buildTXTMatrix('unstructured3D',N,A);
return;
