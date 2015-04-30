%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));

rng(0);
%% Data
load('./../Data/HLC_data.mat');% Xtrain Ytrain Xtest Ytest
[Data] = data_struct(Xtrain, Ytrain, Xtest, Ytest);

%% Parameter - lambda
lam_given = 5 * sqrt(log(Data.p + Data.q) / Data.n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [10 1 0.1 0.01 0.001 1e-4 1e-5 1e-6 1e-7];
    %     lambda_seq = 1;
end

%% CV
kcv = 5; % k-fold CV

%% Optimization method
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
opt_algs = {'PNOPT'};

%% run!
for alg_idx = 1: length(opt_algs)
    alg = opt_algs{alg_idx};
    
    % serach for opt lambda
    opt_para{alg_idx} = ParaOpt(Data, lambda_seq, kcv, alg);  
    opt_lambda = opt_para{alg_idx}.lambda;
    
    % train the model
    opt{alg_idx} = TrainPGM(Data, alg, opt_lambda);
    
    % test the model
    testerr(alg_idx) = PGM_predict(opt{alg_idx}.theta, opt{alg_idx}.alpha1, opt{alg_idx}.beta, opt{alg_idx}.betad, Data.X_te, Data.D_te);
    disp(['Test err : ' num2str(testerr(alg_idx))]);
end



