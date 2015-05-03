%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));

rng(0);
%%
load('./../Data/HLC_data.mat');% Xtrain Ytrain Xtest Ytest
[Data] = data_struct(Xtrain, Ytrain, Xtest, Ytest);


%% para setting
lam_given = 5 * sqrt(log(Data.p + Data.q) / Data.n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [10 1 0.1 0.01 0.001 1e-4 1e-5 1e-6 1e-7];
    %     lambda_seq = 1;
end

kcv = 5; % k-fold CV


%% select optimization algorithms to run
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
opt_algs = {'PNOPT'};


%% run!
opt = TrainPGM_HA(Data, lambda_seq, kcv, opt_algs);


for i = 1 : length(opt_algs)
    testerr(i) = PGM_predict(opt{i}.theta, opt{i}.alpha1, opt{i}.beta, opt{i}.betad, Data.X_te, Data.D_te);
    disp(['Test err : ' num2str(testerr(i))]);
end

