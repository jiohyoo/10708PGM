%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));


%%
load('./../Data/ToyData.mat', 'X', 'Y', 'D', 'p', 'q', 'L', 'n', 'thcts', 'maskDisCts', 'maskDis');
ToyData.p = p;
ToyData.q = q;
ToyData.L = L;

ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.maskDis = maskDis;

X_tr = X;
Y_tr = Y;
D_tr = D;



%% para setting
lam_given = 5 * sqrt(log(ToyData.p + ToyData.q) / n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [lam_given 0.16 0.08 0.04 0.02 0.01 0.005];
end

kcv = 5; % k-fold CV
n_rep = 5;


%% select optimization algorithms to run
opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
%  opt_algs = {'PNOPT'};





%%
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
    
    
    
    %% run!
    opt = TrainPGM(ToyData, lambda_seq, kcv, opt_algs);
    
    for opt_alg_idx = 1: size(opt_algs,2)
        fprintf('---------------------------------------\n');
        [prederr(opt_alg_idx, rep)] = PGM_predict(opt{opt_alg_idx}.theta, opt{opt_alg_idx}.alpha1, opt{opt_alg_idx}.beta, opt{opt_alg_idx}.betad, ToyData.X_te, ToyData.D_te);
    end
end



prederr_opt_alg = mean(prederr,2);
for opt_alg_idx = 1: size(opt_algs,2)
    fprintf('Opt %s (conditional) - lambda: %g\n: prediction error: %g\n', opt_algs{opt_alg_idx}, opt{opt_alg_idx}.lambda, prederr_opt_alg(opt_alg_idx));
    
    if strcmp(opt_algs{opt_alg_idx}, 'PNOPT')
        obj_val(opt_alg_idx) = opt{opt_alg_idx}.out;
    else
        obj_val(opt_alg_idx) = opt{opt_alg_idx}.out.f(end);
    end
end



% save ToyRes




