function main()
clc
make
clc

dic_type 	= 'Rand_sign';
n     		= 512;
redundancy 	= 4;
sigma 		= 0.02;

T 			= round(0.1 * n);
lambda 		= 4 * sigma;

A 		 	= getDictionary(dic_type, n, redundancy);
max_iter 	= 1000;
accuracy 	= 1e-4;

% experiment.clean_signal = b;
% experiment.noise 		  = no;
% experiment.sigma 		  = sigma;
% experiment.supp 		  = supp;
% experiment.x_orig  - exact solution for comparison later

experiment 	= GenerateSyntheticExperiment(A, T, sigma);
y 			= experiment.noisy_signal;

% now we solve min 0.5*||A*x-y||_2^2 + lambda*||x||_1
tic
r = y;

for k=1:1
	x = zeros(size(A,2),1);
	r = y;
	[x1,r,L2L1_CD_iter,atr,fx_trace1,time_trace1,trace_iter1] = solveL2L1_CD_mex(A, x, r, lambda, max_iter, accuracy);
end
L2L1_CD_iter
x1_L2L1_CD = toc


x = zeros(size(A,2),1);
M = A'*A;
c = A'*y;
Mx = zeros(size(A,2),1);

tic
for k=1:1
	[x2,temp,M_CD_iter,fx_trace2,time_trace2,trace_iter2] = solve_M_CD_mex(M, x, Mx, c, lambda, max_iter, accuracy);
end
M_CD_iter
x2_M_CD = toc
fx_trace2 = fx_trace2 + 0.5*y'*y;

ratio = 0.5;
Atr_prev = A' * y;
x = zeros(size(A,2),1);
r = y;

tic
for k=1:1
	[L2L1_MLCD_iter, x3,atr3,r3,fx_trace3,time_trace3,trace_iter3] = solveL2L1_MLCD_mex(A, x, r, lambda, max_iter, accuracy, ratio, Atr_prev);
end
L2L1_MLCD_iter
x3_L2L1_MLCD = toc
trace_iter3




x = zeros(size(A,2),1);
Mx = zeros(size(A,2),1);

tic
for k=1:1
	[M_MLCD_iter,x4,Mx_tmp,fx_trace4,time_trace4,trace_iter4] = solve_M_MLCD_mex(M, x, Mx, c, lambda, max_iter, accuracy, ratio);
end
M_MLCD_iter
x4_M_MLCD = toc
fx_trace4 = fx_trace4 + 0.5*y'*y;

opt = min([min(fx_trace1(1:trace_iter1)),...
           min(fx_trace3(1:trace_iter3)),...
           min(fx_trace2(1:trace_iter2)),...
           min(fx_trace4(1:trace_iter4))]); 


figure;
semilogy(time_trace1(1:trace_iter1),fx_trace1(1:trace_iter1) - opt,'b-*','LineWidth',1.6);
hold on;
semilogy(time_trace3(1:trace_iter3),fx_trace3(1:trace_iter3) - opt,'r-*','LineWidth',1.6);
semilogy(time_trace2(1:trace_iter2),fx_trace2(1:trace_iter2) - opt,'k-*','LineWidth',1.6)
semilogy(time_trace4(1:trace_iter4),fx_trace3(1:trace_iter4) - opt,'m-*','LineWidth',1.6)
legend('A^T*A_{CD}','A^T*A_{MLCD}','M_{CD}','M_{MLCD}');
title('f(x)-f_{opt} history: time','FontSize',14);
xlabel('seconds','FontSize',14)

norm_x1x2 = norm(x1-x2)
norm_x3x1 = norm(x3-x1)
norm_x4x1 = norm(x4-x1)
norm_x4x3 = norm(x4-x3)

dic_type
