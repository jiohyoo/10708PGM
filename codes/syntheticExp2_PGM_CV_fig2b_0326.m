%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));



%% load data
load ToyData X Y D p q L n thcts maskDisCts maskDis
% load ToyData_v2 ToyData
ToyData.p = p;
ToyData.q = q;
ToyData.L = L;
% ToyData.n = n;
ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.maskDis = maskDis;
% ToyData.X_tr = X;
% ToyData.Y_tr = Y;
% ToyData.D_tr = D;
% ToyData.X_te = [];
% ToyData.Y_te = [];
% ToyData.D_te = [];

X_tr = X;
Y_tr = Y;
D_tr = D;

% n_interval = 250;
% n_step = 1;%ToyData.n / n_interval;
n_rep = 10;

%% para setting
lam_given = 5 * sqrt(log(ToyData.p + ToyData.q) / n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [1 0.64 0.32 lam_given 0.16 0.08 0.04 0.02];
end



for rep = 1: n_rep
    disp(num2str(rep));
    
    [TRAIN, TEST] = crossvalind('LeaveMOut', n ,400);
    X_tr_CV = X_tr(find(TRAIN), :);
    D_tr_CV = D_tr(find(TRAIN), :);
    Y_tr_CV = Y_tr(find(TRAIN),:);
    X_te_CV  = X_tr(find(TEST), :);
    Y_te_CV  = Y_tr(find(TEST), :);
    D_te_CV  = D_tr(find(TEST), :);
    
    
    %     ToyData.n = n;
    ToyData.X_tr = X_tr_CV;
    ToyData.Y_tr = Y_tr_CV;
    ToyData.D_tr = D_tr_CV;
    ToyData.X_te = X_te_CV;
    ToyData.Y_te = Y_te_CV;
    ToyData.D_te = D_te_CV;
    
    ToyData.n_tr = 1600;
    
    
    opt = syntheticExp2_test_CV_fn(ToyData, lambda_seq, 5);
    
    %% evaluate prediction error
    fprintf('---------------------------------------\n');
    [prederr(rep)] = test_predict(opt.theta, opt.alpha1, opt.beta, opt.betad, ToyData.Y_te, ToyData.D_te);
    fprintf('Original_model - lambda: %g\n, prediction error: %g\n', opt.lambda, prederr(rep));
    
    
end
