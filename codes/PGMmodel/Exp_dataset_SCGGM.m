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
% load('./../HLC_data_Toy.mat');% Xtrain Ytrain Xtest Ytest
%  load('/HLC_data_Toy.mat');% Xtrain Ytrain Xtest Ytest
Dir = ['C:\Users\HyunAh\Desktop\10708PGM\codes\SCGGM_code\data\'];
xtrain = load([Dir 'xtrain.txt']);
ytrain = load([Dir 'ytrain.txt']);

xtest = load([Dir 'xtest.txt']);
ytest = load([Dir 'ytest.txt']);
% [Data] = data_struct(Xtrain_toy, Ytrain_toy, Xtest_toy, Ytest_toy);


Data.n = size(xtrain, 1);
Data.n_tr = Data.n;
Data.n_te = size(xtest, 1);


Data.Y_tr = xtrain;
Data.X_tr = ytrain;
%
Data.p = size(ytrain,2);
Data.q = size(xtrain,2);
Data.L = 3*ones(Data.q ,1);

Dtr=[];
for j=1:Data.q 
    Dj=zeros(Data.n_tr,Data.L(j));
    for i=1:Data.n_tr
        Dj(i,xtrain(i,j)+1)=1;
    end
    Dtr=[Dtr Dj];
end
Data.D_tr = Dtr;
Data.Y_te = xtest;
Data.X_te = ytest;
% 

Dte=[];
for j=1:Data.q 
    Dj=zeros(Data.n_te,Data.L(j));
    for i=1:Data.n_te
        Dj(i,xtest(i,j)+1)=1;
    end
    Dte=[Dte Dj];
end
Data.D_te = Dte;


%% Parameter - lambda
lam_given = 5 * sqrt(log(Data.p + Data.q) / Data.n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [1e-3 1e-5 1e-7];
    %     lambda_seq = 1;
end

%% CV
kcv = 3; % k-fold CV

%% Optimization method
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
opt_algs = {'PNOPT'};

%% run!
for alg_idx = 1: length(opt_algs)
    alg = opt_algs{alg_idx};
    
    % serach for opt lambda
    opt_para{alg_idx} = ParaOpt(Data, lambda_seq, kcv, alg);  
    opt_lambda = opt_para{alg_idx}.lambda;
    
    disp(['opt lambda: ' num2str(opt_lambda)]);
    
    % train the model
    opt{alg_idx} = TrainPGM(Data, alg, opt_lambda);
    
    % test the model
    testerr(alg_idx) = PGM_predict(opt{alg_idx}.theta, opt{alg_idx}.alpha1, opt{alg_idx}.beta, opt{alg_idx}.betad, Data.X_te, Data.D_te);
    disp(['Test err : ' num2str(testerr(alg_idx))]);
end



