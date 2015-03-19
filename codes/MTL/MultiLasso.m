% Multi task Lasso with L1/L2 regularization
% TODO List:
% 1) LOAD DATA- Fill in X, Y
% 2) Set Train/Test set ratio
% 3) Set cv_fold (== 5 in Sohn_2012)
% 4) Prediction error measure (MSE?) 


clear; clc;

addpath(genpath('.'));


X = load_data.X;
Y = load_data.Y;

% preprocessing data
for t = 1: length(X)
    X{t} = zscore(X{t});                  % normalization
    X{t} = [X{t} ones(size(X{t}, 1), 1)]; % add bias. 
end

% split data into training and testing.
training_percent = 0.3;
[X_tr, Y_tr, X_te, Y_te] = mtSplitPerc(X, Y, training_percent);

% the function used for evaluation.
eval_func_str = 'eval_MTL_mse';
higher_better = false;  % mse is lower the better.

% cross validation fold
cv_fold = 5;

% optimization options
opts = [];
opts.maxIter = 100;

% model parameter range
param_range = [0.001 0.01 0.1 1 10 100 1000 10000];

fprintf('Perform model selection via cross validation: \n')
[ best_param perform_mat] = CrossValidation1Param...
    ( X_tr, Y_tr, 'Least_L21', opts, param_range, cv_fold, eval_func_str, higher_better);

%disp(perform_mat) % show the performance for each parameter.

% build model using the optimal parameter 
W = Least_L21(X_te, Y_te, best_param, opts);

% show final performance
final_performance = eval_MTL_mse(Y_te, X_te, W);
fprintf('Performance on test data: %.4f\n', final_performance);
