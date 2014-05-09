function experiment = GenerateSyntheticExperiment(A, T, sigma)

% T - number of atoms in the support
% We assume \sigma_x = 1

[N,M] = size(A);
x = zeros(M,1);

supp = randperm(M);
supp = supp(1:T);

x(supp) = sign(randn(T,1));
b = A*x;
b = b/max(abs(b));

no = sigma*randn(N,1);
y = b + no;
experiment.m = size(A,2);
%experiment.A = A;
experiment.clean_signal = b;
experiment.noise = no;
experiment.sigma = sigma;
experiment.noisy_signal = y;
experiment.supp = supp;
experiment.x_orig = x;
experiment.n = N;
return;
