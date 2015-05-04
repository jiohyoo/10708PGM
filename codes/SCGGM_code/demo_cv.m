%--------------------------------------------------------------------------
% Sample run of the sparse CGGM algorithm with cross-validation
%--------------------------------------------------------------------------

%% sCGGM demo with cross-validation 
% specify the search grid of regularization parameters
lambda1_seq = [0.16 0.08 0.04 0.02]; 
lambda2_seq = [0.16 0.08 0.04 0.02]; 

% performs kcv-fold cross validation, kcv must be >= 3
kcv = 5;  

% loading traing data and test data 
xtrain = load('./data/xtrain.txt');
ytrain = load('./data/ytrain.txt');

xtest = load('./data/xtest.txt');
ytest = load('./data/ytest.txt');

fprintf('sCGGM demo with %d-fold cross-validation...\n', kcv);

option.verbose = true; 
option.maxiter = 500; 

opt = scggm_cv( xtrain, ytrain, kcv, lambda1_seq, lambda2_seq, option);

% compute prediction errors
[~, e] = scggm_predict(opt.Theta, opt.intercept,  xtest, ytest);
fprintf('sCGGM demo  completed, test set prediction error: %g\n', e); 

% perform inference
[Beta] = scggm_indirect_SNP_overall(opt.Theta); 

% decomposition of overall indirect SNP perturbations
% passed by the k-th gene
k = 2;  % k = 1 ... 30 
[Beta_k] = scggm_indirect_SNP_decompose(opt.Theta, k);

% decomposition of gene-expression covariance
Cov = scggm_cov_decompose(opt.Theta, xtrain, ytrain);

if ~exist('results/demo_cv', 'dir')
	mkdir('results/demo_cv'); 
end

dlmwrite('./results/demo_cv/optimal_Theta_xy.txt', opt.Theta.xy, '\t');
dlmwrite('./results/demo_cv/optimal_Theta_yy.txt', opt.Theta.yy, '\t'); 
dlmwrite('./results/demo_cv/optimal_intercept.txt', opt.intercept, '\t');
dlmwrite('./results/demo_cv/optimal_lambdas.txt', opt.lambdas, '\t');
dlmwrite('./results/demo_cv/Beta.txt', Beta, '\t');
dlmwrite(['./results/demo_cv/Beta_', num2str(k), '.txt'], Beta_k, '\t');
dlmwrite('./results/demo_cv/Cov_Overall.txt', Cov.Overall, '\t'); 
dlmwrite('./results/demo_cv/Cov_Network_Induced.txt', Cov.Network_Induced, '\t');
dlmwrite('./results/demo_cv/Cov_SNP_Induced.txt', Cov.SNP_Induced, '\t');


figure;imagesc(opt.Theta.xy);
figure;imagesc(opt.Theta.yy);