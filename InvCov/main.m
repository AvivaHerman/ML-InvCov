function [] = main()
close all
clc

n = 4096; % matrix size: nxn.
% k = 1*n; 
k = 512;%round(20*sqrt(n));  % number of samples in X.

rng(1); % set random seed

fprintf('Generating Graph...');
[A_true,X] = generateExperiment(n,k);
S = ((1/k)*X)*X';
clear X;

fprintf('Done\n');


% load('ER_692.mat');

lambda = 0.2;
params.epsilon = 1e-2;
%params.epsilon_threshold = 1e-7;
make

mode = 'trace';
msg_urgency = 2;
max_iter = 20;



[~ , ~ , opt_res1, time_res1, num_iter1, supp_res1, numActive_res1] = QUIC(mode, S, lambda, params.epsilon, msg_urgency, max_iter);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~QUIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
pause(0.5);
[A , ~, opt_res2, time_res2, num_iter2, supp_res2, numActive_res2, num_vcycleIter2] = ML_QUIC(mode, S, lambda, params.epsilon, 2, max_iter);
% 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ML_QUIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');


% figure;
% subplot(1,2,1);
% spy(A_true);
% title('A-true')
% subplot(1,2,2);
% spy(A);
% title('A')


% A = A.*(abs(A)>1e-2*max(max(abs(A))));
% figure;
% subplot(1,2,1);
% spy(A_true);
% title('A-true')
% subplot(1,2,2);
% spy(A);
% title('A')


% 

num_iter2 = max(find(opt_res2>0));
num_vcycles = max(find(num_vcycleIter2>0));

opt = min(opt_res1(num_iter1),opt_res2(num_iter2));

% figure;
% semilogy(1:num_iter1,opt_res1(1:num_iter1) - opt,'b-*','LineWidth',1.6);
% hold on
% semilogy(1:num_iter2,opt_res2(1:num_iter2) - opt,'r-','LineWidth',1.6);
% semilogy(num_vcycleIter2(1:num_vcycles),opt_res2(num_vcycleIter2(1:num_vcycles)) - opt,'kp','LineWidth',1.6);
% semilogy(1:num_iter3,opt_res3(1:num_iter3) - opt,'m-*','LineWidth',1.6);
% legend('QUIC','MLQUIC','MLQUIC-Vcycle end','QUIC-original');
% title('f(x)-f_{opt} history: iterations','FontSize',14);
% xlabel('num relaxations','FontSize',14)
% 
% 
figure;
semilogy(time_res1(1:num_iter1),opt_res1(1:num_iter1) - opt,'b-*','LineWidth',1.6);
hold on
semilogy(time_res2(1:num_iter2),opt_res2(1:num_iter2) - opt,'r-','LineWidth',1.6);
semilogy(time_res2(num_vcycleIter2(1:num_vcycles)),opt_res2(num_vcycleIter2(1:num_vcycles)) - opt,'kp');
%semilogy(time_res3(1:num_iter3),opt_res3(1:num_iter3) - opt,'m-*');
legend('QUIC','MLQUIC','MLQUIC-Vcycle end','QUIC-original');
title('f(x)-f_{opt} history: time','FontSize',14);
xlabel('seconds','FontSize',14)

figure;
plot(time_res1(1:num_iter1),supp_res1(1:num_iter1),'b','LineWidth',1.6);
hold on
plot(time_res2(1:num_iter2),supp_res2(1:num_iter2),'r','LineWidth',1.6);
plot(time_res2(num_vcycleIter2(1:num_vcycles)),supp_res2(num_vcycleIter2(1:num_vcycles)),'kp');
plot(time_res1(1:num_iter1),numActive_res1(1:num_iter1),':b','LineWidth',1.6);
plot(time_res2(1:num_iter2),numActive_res2(1:num_iter2),':r','LineWidth',1.6);
plot(time_res2(num_vcycleIter2(1:num_vcycles)),numActive_res2(num_vcycleIter2(1:num_vcycles)),'kp');
legend('QUIC-Support','MLQUIC-Support','MLQUIC-Vcycle end','QUIC-ActiveSet','MLQUIC-ActiveSet');
title('Support and Free Set Sizes History: time','FontSize',14);
xlabel('seconds','FontSize',14)

figure;
plot(2:num_iter1,supp_res1(2:num_iter1),'b','LineWidth',1.6);
hold on
plot(2:num_iter2,supp_res2(2:num_iter2),'r','LineWidth',1.6);
plot(num_vcycleIter2(2:num_vcycles),supp_res2(num_vcycleIter2(2:num_vcycles)),'kp','LineWidth',1.6);
plot(2:num_iter1,numActive_res1(2:num_iter1),':b','LineWidth',1.6);
plot(2:num_iter2,numActive_res2(2:num_iter2),':r','LineWidth',1.6);
plot(num_vcycleIter2(2:num_vcycles),numActive_res2(num_vcycleIter2(2:num_vcycles)),'kp');
legend('QUIC-Support','MLQUIC-Support','MLQUIC-Vcycle end', 'QUIC-ActiveSet','MLQUIC-ActiveSet');
title('Support and Active Set Sizes History: iterations','FontSize',14);
xlabel('num relaxations','FontSize',14)





