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
ToyData.n = n;
ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.maskDis = maskDis;
ToyData.X = X;
ToyData.Y = Y;
ToyData.D = D;
ToyData.X_te = [];
ToyData.Y_te = [];
ToyData.D_te = [];

n_interval = 250;
n_step = ToyData.n / n_interval;
n_rep = 10;


for ii = 1: n_step
    
    n_sample = n_interval * ii;
    ToyData.n_tr = n_sample;
    
    for rep = 1: n_rep
        
        rand_perm = randperm(ToyData.n);
        ToyData.X_tr = ToyData.X(rand_perm(1:n_sample), :);
        ToyData.Y_tr = ToyData.Y(rand_perm(1:n_sample), :);
        ToyData.D_tr = ToyData.D(rand_perm(1:n_sample), :);
        
        
        
        
        %% para setting
        lam_given = 5 * sqrt(log(ToyData.p + ToyData.q) / ToyData.n_tr);
        use_given_lam = 0;
        
        if use_given_lam
            lambda_seq = lam_given;
        else
            lambda_seq = [1 0.64 0.32 lam_given 0.16 0.08 0.04 0.02];
        end
        
        
        
        %% run!
        kcv = 5; % k-fold CV
        opt = syntheticExp2_PGM_CV_fn(ToyData, lambda_seq, kcv);
        
        
        
        %% plot graph!
%         opt.theta' % cts-dis 10 by 20
%         ToyData.maskDisCts % cts-dis truth
        theta_edges = (opt.theta' > 0); % binary
        theta_edges2 = zeros(size(theta_edges,1), size(theta_edges,2)/2);
        for jj = 1: size(theta_edges,2)/2
           theta_edges2(:,jj) = sum(theta_edges(:,(jj-1)*2+1:(jj-1)*2+2),2); % 2 states - OR op
        end
        theta_edges2 = theta_edges2 > 0;
%         n_corr_rec = sum(sum(theta_edges2 .* ToyData.maskDisCts));
%         perc_corr_rec(ii, rep) = n_corr_rec / sum(sum(ToyData.maskDisCts));
        TP(ii,rep) = sum(sum(theta_edges2 .* ToyData.maskDisCts));
        FP(ii,rep) = sum(sum(theta_edges2 .* (ToyData.maskDisCts * (-1) + 1)));
        TN(ii,rep) = sum(sum((theta_edges2 * (-1) + 1) .* (ToyData.maskDisCts * (-1) + 1)));
        FN(ii,rep) = sum(sum((theta_edges2 * (-1) + 1) .* ToyData.maskDisCts));
        
    end
end

%% plot!
% perc_corr_rec_avg = mean(perc_corr_rec,2);
% figure;
% plot(perc_corr_rec_avg);
TP_avg = mean(TP,2);
FP_avg = mean(FP,2);
TN_avg = mean(TN,2);
FN_avg = mean(FN,2);
Accuracy = (TP_avg + TN_avg) ./ (TP_avg + FP_avg + TN_avg + FN_avg);
figure;
plot([0:250:2000], [0; Accuracy]);

save plot2b TP FP TN FN TP_avg FP_avg TN_avg FN_avg Accuracy
 

%% evaluate prediction error
% fprintf('---------------------------------------\n');
% [prederr] = PGM_predict(opt.theta, opt.alpha1, opt.beta, opt.betad, ToyData.Y_te, ToyData.D_te);
% fprintf('PGM_model - lambda: %g\n, prediction error: %g\n', opt.lambda, prederr);



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



 