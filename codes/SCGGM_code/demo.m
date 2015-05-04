%--------------------------------------------------------------------------
% Sample run of the sparse CGGM algorithm without cross-validation
%--------------------------------------------------------------------------
% specify regularization parameters
lambda_1 = 0.1; 
lambda_2 = 0.1; 

% loading dataset 
xtrain = load('./data/xtrain.txt');
ytrain = load('./data/ytrain.txt');
fprintf('sCGGM demo...\nJ = %d, K = %d, sample size = %d\nRegularization parameters: lambda_1 = %g, lambda_2 = %g\n', size(xtrain,2), size(ytrain,2), size(xtrain,1), lambda_1, lambda_2);

% run SCGGM
opt = scggm(xtrain, ytrain, lambda_1, lambda_2); 

% overall indirect SNP perturbations
Beta = scggm_indirect_SNP_overall(opt.Theta); 

% decomposition of overall indirect SNP perturbations
% passed by the k-th gene
k = 2; 
Beta_k = scggm_indirect_SNP_decompose(opt.Theta, k);

% decomposition of gene-expression covariance
Cov = scggm_cov_decompose(opt.Theta, xtrain, ytrain);

if ~exist('results/demo', 'dir')
	mkdir('results/demo'); 
end
dlmwrite('./results/demo/optimal_Theta_xy.txt', opt.Theta.xy, '\t');
dlmwrite('./results/demo/optimal_Theta_yy.txt', opt.Theta.yy, '\t'); 
dlmwrite('./results/demo/optimal_intercept.txt', opt.intercept, '\t');
dlmwrite('./results/demo/Beta.txt', Beta, '\t');
dlmwrite(['./results/demo/Beta_', num2str(k), '.txt'], Beta_k, '\t');
dlmwrite('./results/demo/Cov_Overall.txt', Cov.Overall, '\t'); 
dlmwrite('./results/demo/Cov_Network_Induced.txt', Cov.Network_Induced, '\t');
dlmwrite('./results/demo/Cov_SNP_Induced.txt', Cov.SNP_Induced, '\t');



[~, e] = scggm_predict(opt.Theta, opt.intercept, X_te, Y_te);
figure;imagesc(opt.Theta.xy);
figure;imagesc(opt.Theta.yy);
