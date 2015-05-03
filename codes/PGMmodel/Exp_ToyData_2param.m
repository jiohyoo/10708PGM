%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));

%% Data
%load('./../Data/ToyData4states.mat');% Xtrain Ytrain Xtest Ytest
%[Data] = data_struct(Xtrain, Ytrain, Xtest, Ytest);

%% Data
load('./../Data/ToyData4states.mat', 'X', 'Y', 'D', 'p', 'q', 'L', 'n', 'thcts', 'maskDisCts', 'maskDis', 'thdiscts');
ToyData.p = p;
ToyData.q = q;
ToyData.L = L;
ToyData.n_tr = n;
ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.thdiscts = thdiscts;
ToyData.maskDis = maskDis;
ToyData.X_tr = X;
ToyData.Y_tr = Y;
ToyData.D_tr = D;
ToyData.X_te = X;
ToyData.Y_te = Y;
ToyData.D_te = D;
Data = ToyData

%% Parameter - lambda1, lambda2
lambda1_seq = [1e-2 1e-4 1e-6];
lambda2_seq = [1e-2 1e-4 1e-6];

%% CV
kcv = 3; % k-fold CV

%% Optimization method
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
opt_algs = {'PNOPT'};

%% run!
for alg_idx = 1: length(opt_algs)
    alg = opt_algs{alg_idx};
    
    % serach for opt lambda
    opt_para{alg_idx} = ParaOpt_2param(Data, lambda1_seq, lambda2_seq, kcv, alg);  
    opt_lambda1 = opt_para{alg_idx}.lambda1;
    opt_lambda2 = opt_para{alg_idx}.lambda2;
    
    % train the model
    opt{alg_idx} = TrainPGM_2param(Data, alg, opt_lambda1, opt_lambda2);
    
    % test the model
    testerr(alg_idx) = PGM_predict(opt{alg_idx}.theta, opt{alg_idx}.alpha1, opt{alg_idx}.beta, opt{alg_idx}.betad, Data.X_te, Data.D_te);
    disp(['Test err : ' num2str(testerr(alg_idx))]);
end

