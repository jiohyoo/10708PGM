%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));

%% Data
% load('./../Data/ToyData.mat', 'X', 'Y', 'D', 'p', 'q', 'L', 'n', 'thcts', 'maskDisCts', 'maskDis');
% load ToyData2_n1000.mat
load ToyData2_q300_n1000.mat
ToyData.p = p;
ToyData.q = q;
ToyData.L = L;

ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.maskDis = maskDis;

X_tr = X;
Y_tr = Y;
D_tr = D;

%% Parameter - lambda
lam_given = 5 * sqrt(log(ToyData.p + ToyData.q) / n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [lam_given 0.1 0.01];
end

%% CV
kcv = 5; % k-fold CV
n_rep = 1;

%% Optimization method
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
 opt_algs = {'PNOPT'};

%% run
for alg_idx = 1: length(opt_algs)
    alg = opt_algs{alg_idx};

    for rep = 1: n_rep
        disp(num2str(rep));
        
        [TRAIN, TEST] = crossvalind('LeaveMOut', n , n/kcv);
        X_tr_CV = X_tr(find(TRAIN), :);
        D_tr_CV = D_tr(find(TRAIN), :);
        Y_tr_CV = Y_tr(find(TRAIN),:);
        X_te_CV  = X_tr(find(TEST), :);
        Y_te_CV  = Y_tr(find(TEST), :);
        D_te_CV  = D_tr(find(TEST), :);
        
        ToyData.X_tr = X_tr_CV;
        ToyData.Y_tr = Y_tr_CV;
        ToyData.D_tr = D_tr_CV;
        ToyData.X_te = X_te_CV;
        ToyData.Y_te = Y_te_CV;
        ToyData.D_te = D_te_CV;
        
        ToyData.n_tr = n - n / kcv;
        
 
        % serach for opt lambda
        opt_para = ParaOpt(ToyData, lambda_seq, kcv, alg);
        opt_lambda = opt_para.lambda;
%         
        % train the model
%         opt_lambda = 0.001;
        opt = TrainPGM(ToyData, alg, opt_lambda);
        
        % test the model
        fprintf('---------------------------------------\n');
        testerr(alg_idx, rep) = PGM_predict(opt.theta, opt.alpha1, opt.beta, opt.betad, ToyData.X_te, ToyData.D_te);
        disp(['Alg-' alg ' # rep-' num2str(rep) ' - Test err : ' num2str(testerr(alg_idx, rep))]);
    end
end



