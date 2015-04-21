%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));


%%
load HLC_data.mat Xtrain Ytrain Xtest Ytest
Data.p = size(Ytrain,2);
Data.q = size(Xtrain,2);
Data.n = size(Xtrain,1);
Data.L = 3*ones(size(Xtrain,2),1);

Dtr = zeros(Data.n, Data.L(1) * Data.p);
for i = 1: Data.n
    for j = 1: Data.q
        Dtr(i, (j-1)*Data.L(1)+Xtrain(i,j)+1) = 1;
    end
end

Dte = zeros(size(Xtest,1), Data.L(1) * size(Xtest,2));
for i = 1: size(Xtest,1)
    for j = 1: size(Xtest,2)
        Dte(i, (j-1)*Data.L(1)+Xtest(i,j)+1) = 1;
    end
end

Data.Y_tr = Xtrain;
Data.X_tr = Ytrain;
Data.D_tr = Dtr;
Data.Y_te = Xtest;
Data.X_te = Ytest;
Data.D_te = Dte;

Data.n_tr = size(Ytrain,1);





%% para setting
lam_given = 5 * sqrt(log(Data.p + Data.q) / Data.n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [lam_given 0.16 0.08 0.04 0.02 0.01 0.005];
end

kcv = 5; % k-fold CV



%% select optimization algorithms to run
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
opt_algs = {'PNOPT'};



%%
rep = 1;

    %% run!
    opt = syntheticExp2_PGM_CV_fn_opts(Data, lambda_seq, kcv, opt_algs);
    
    % original algorithm - learns joint distribution
    % % %     opt_ori = syntheticExp2_PGM_CV_fn_original(ToyData, lambda_seq, kcv);
    
    
    
    %     %% evaluate prediction error
    %     fprintf('---------------------------------------\n');
    %     [prederr(rep)] = PGM_predict(opt.theta, opt.alpha1, opt.beta, opt.betad, ToyData.Y_te, ToyData.D_te);
    %     fprintf('PGM_model (conditional) - lambda: %g\n, prediction error: %g\n', opt.lambda, prederr(rep));
    % % %
    % % %     fprintf('---------------------------------------\n');
    % % %     [prederr_ori(rep)] = PGM_predict(opt_ori.theta, opt_ori.alpha1, opt_ori.beta, opt_ori.betad, ToyData.Y_te, ToyData.D_te);
    % % %     fprintf('PGM_model (joint) - lambda: %g\n, prediction error: %g\n', opt_ori.lambda, prederr_ori(rep));
    for opt_alg_idx = 1: size(opt_algs,2)
        fprintf('---------------------------------------\n');
        [prederr(opt_alg_idx, rep)] = PGM_predict(opt{opt_alg_idx}.theta, opt{opt_alg_idx}.alpha1, opt{opt_alg_idx}.beta, opt{opt_alg_idx}.betad, Data.X_te, Data.D_te);
%         fprintf('PGM_model (conditional) - lambda: %g\n, prediction error: %g\n', opt{opt_alg_idx}.lambda, prederr(opt_alg_idx,rep));
    end
% end



prederr_opt_alg = mean(prederr,2);
for opt_alg_idx = 1: size(opt_algs,2)
    fprintf('Opt %s (conditional) - lambda: %g\n: prediction error: %g\n', opt_algs{opt_alg_idx}, opt{opt_alg_idx}.lambda, prederr_opt_alg(opt_alg_idx));
    
    if strcmp(opt_algs{opt_alg_idx}, 'PNOPT')
        obj_val(opt_alg_idx) = opt{opt_alg_idx}.out;
    else
         obj_val(opt_alg_idx) = opt{opt_alg_idx}.out.f(end);
    end
end



% save ResOut




