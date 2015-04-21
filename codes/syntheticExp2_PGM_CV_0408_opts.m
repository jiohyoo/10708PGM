%% Demo of the synthetic experiments.
%% Requires UGM at http://www.di.ens.fr/~mschmidt/Software/UGM_2009.zip
%% Requires TFOCS at http://tfocs.stanford.edu
clear all
close all

addpath(genpath('./UGM_2009'));
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));


%%
load ToyData X Y D p q L n thcts maskDisCts maskDis
% load ToyData_v2 ToyData
ToyData.p = p;
ToyData.q = q;
ToyData.L = L;

ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.maskDis = maskDis;

X_tr = X;
Y_tr = Y;
D_tr = D;

% n_interval = 250;
% n_step = 1;%ToyData.n / n_interval;
n_rep = 1;



%% para setting
lam_given = 5 * sqrt(log(ToyData.p + ToyData.q) / n);
use_given_lam = 0;

if use_given_lam
    lambda_seq = lam_given;
else
    lambda_seq = [lam_given 0.16 0.08 0.04 0.02 0.01 0.005];
%     lambda_seq = [0.01 0.005 0.0025];
end

kcv = 5; % k-fold CV



%% select optimization algorithms to run
opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
%  opt_algs = {'PNOPT'};





%%
for rep = 1: n_rep
    disp(num2str(rep));
    
    [TRAIN, TEST] = crossvalind('LeaveMOut', n ,400);
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
    
    ToyData.n_tr = 1600;
    
    
    
    %% run!
    opt = syntheticExp2_PGM_CV_fn_opts(ToyData, lambda_seq, kcv, opt_algs);
    
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
        [prederr(opt_alg_idx, rep)] = PGM_predict(opt{opt_alg_idx}.theta, opt{opt_alg_idx}.alpha1, opt{opt_alg_idx}.beta, opt{opt_alg_idx}.betad, ToyData.X_te, ToyData.D_te);
%         fprintf('PGM_model (conditional) - lambda: %g\n, prediction error: %g\n', opt{opt_alg_idx}.lambda, prederr(opt_alg_idx,rep));
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



% save ResOut_0420_2




%% Plot parameters
% if use_given_lam
%     close all;
%     data_name = 'fixed-lam -cts truth';
%     figure(1); imagesc(triu(ToyData.thcts-diag(diag(ToyData.thcts)))); title(data_name); colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     data_name = 'fixed-lam -cts recover';
%     figure(2); imagesc(-opt.beta); title(data_name); colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     data_name = 'fixed-lam -discrete truth';
%     figure(3); imagesc(ToyData.maskDis); title(data_name);colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     % figure; imagesc(phi); title('dis recover');colorbar;
%     data_name = 'fixed-lam -cts-dis truth';
%     figure; imagesc(ToyData.maskDisCts); title(data_name);colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     data_name = 'fixed-lam -cts-dis recover';
%     figure; imagesc(opt.theta'); title(data_name);colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     
% else
%     close all;
%     data_name = 'cts truth';
%     figure(1); imagesc(triu(ToyData.thcts-diag(diag(ToyData.thcts)))); title(data_name); colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     data_name = 'cts recover';
%     figure(2); imagesc(-opt.beta); title(data_name); colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     data_name = 'discrete truth';
%     figure(3); imagesc(ToyData.maskDis); title(data_name);colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     % figure; imagesc(phi); title('dis recover');colorbar;
%     data_name = 'cts-dis truth';
%     figure; imagesc(ToyData.maskDisCts); title(data_name);colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
%     data_name = 'cts-dis recover';
%     figure; imagesc(opt.theta'); title(data_name);colorbar;
%     saveas(gcf,sprintf('./PLOTS/%s_lam_%.2f_pred_%.4f_res.png', data_name, opt.lambda, prederr));
% end



 