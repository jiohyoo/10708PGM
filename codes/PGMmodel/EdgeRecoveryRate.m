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


kcv = 5; % k-fold CV



%% select optimization algorithms to run
% opt_algs = {'AT', 'GRA','LLM','N07','N83','TS','PNOPT'};
opt_algs = {'PNOPT'};


%%
for ii = 1: n_step
    
    n_sample = n_interval * ii;
    ToyData.n_tr = n_sample;
    
    for rep = 1: n_rep
        disp(num2str(rep));
        
        
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
            lambda_seq = [lam_given 0.16 0.08 0.04 0.02 0.01 0.005];
        end
        
        
        %% run!
        opt = TrainPGM(ToyData, lambda_seq, kcv, opt_algs);
   
        theta_edges = (opt{1}.theta' > 0); % binary
        theta_edges2 = zeros(size(theta_edges,1), size(theta_edges,2)/2);
        for jj = 1: size(theta_edges,2)/2
            theta_edges2(:,jj) = sum(theta_edges(:,(jj-1)*2+1:(jj-1)*2+2),2); % 2 states - OR op
        end
        theta_edges2 = theta_edges2 > 0;

        TP(ii,rep) = sum(sum(theta_edges2 .* ToyData.maskDisCts));
        FP(ii,rep) = sum(sum(theta_edges2 .* (ToyData.maskDisCts * (-1) + 1)));
        TN(ii,rep) = sum(sum((theta_edges2 * (-1) + 1) .* (ToyData.maskDisCts * (-1) + 1)));
        FN(ii,rep) = sum(sum((theta_edges2 * (-1) + 1) .* ToyData.maskDisCts));
        
    end
    
end



%% plot!
FP_avg = mean(FP,2);
TN_avg = mean(TN,2);
FN_avg = mean(FN,2);
Accuracy = (TP_avg + TN_avg) ./ (TP_avg + FP_avg + TN_avg + FN_avg);
figure;
plot([0:250:2000], [0; Accuracy]);

% save EdgeRecoveryRes TP FP TN FN TP_avg FP_avg TN_avg FN_avg Accuracy


