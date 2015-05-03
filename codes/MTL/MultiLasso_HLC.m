% Multi task Lasso with L1/L2 regularization
% TODO List:
% 1) LOAD DATA- Fill in X, Y
% 2) Set Train/Test set ratio
% 3) Set cv_fold (== 5 in Sohn_2012)
% 4) Prediction error measure (MSE?) 


clear; clc;

addpath(genpath('.'));

% loading data
hlc = open('../Data/HLC_data.mat')
hlcX = [hlc.Xtrain_raw ; hlc.Xtest_raw];
hlcY = [hlc.Ytrain ; hlc.Ytest];
[n, q] = size(hlcX);
[~, p] = size(hlcY);

n_iter = 1;
performances = zeros(n_iter, 1);

% preprocessing data
for t = 1:p
    %X{t} = zscore(hlcX);                  % normalization
    X{t} = hlcX;                  % normalization
    %X{t} = [X{t} ones(size(X{t}, 1), 1)]; % add bias. 
	Y{t} = hlcY(:, t);
end

% split data into training and testing.
for iter=1:n_iter
training_percent = 0.8;
%n_tr = floor(n * training_percent);
n_tr = 143;
n_te = n - n_tr;
rand_perm = randperm(n);
for t = 1:p
%	X_tr{t} = X{t}(rand_perm(1:n_tr), :);
%	X_te{t} = X{t}(rand_perm(n_tr+1:end), :);

%	Y_tr{t} = Y{t}(rand_perm(1:n_tr));
%	Y_te{t} = Y{t}(rand_perm(n_tr+1:end));
	X_tr{t} = X{t}(1:143, :);
	X_te{t} = X{t}(144:end, :);

	Y_tr{t} = Y{t}(1:143);
	Y_te{t} = Y{t}(144:end);
end


% the function used for evaluation.
eval_func_str = 'eval_mse';
higher_better = false;  % mse is lower the better.

% cross validation fold
cv_fold = 5;

% optimization options
opts = [];
opts.maxIter = 100;

% model parameter range
%param_range = [0.001 0.01 0.1 1 10 100 1000 10000];
param_range = [0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56]

fprintf('Perform model selection via cross validation: \n')
[ best_param perform_mat] = CrossValidation1Param...
    ( X_tr, Y_tr, 'Least_L21', opts, param_range, cv_fold, eval_func_str, higher_better);

%disp(perform_mat) % show the performance for each parameter.

% build model using the optimal parameter 
W = Least_L21(X_te, Y_te, best_param, opts);

% show final performance
final_performance = eval_mse(Y_te, X_te, W);
fprintf('Performance on test data: %.4f\n', final_performance);
performances(iter) = final_performance;
end
