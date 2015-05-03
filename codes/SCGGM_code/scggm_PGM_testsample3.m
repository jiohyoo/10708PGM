%--------------------------------------------------------------------------
% Running SCGGM with ToyData  
% 1) 5-fold CV, 2) prediction error
%--------------------------------------------------------------------------

%% sCGGM demo with cross-validation 
% specify the search grid of regularization parameters
lambda1_seq = [0.04 0.02 0.01 0.005]; 
lambda2_seq = [0.04 0.02 0.01 0.005]; 

% performs kcv-fold cross validation, kcv must be >= 3
kcv = 5;  

% loading data
% data = load('./../ToyData.mat');
data = load('./../ToyData3.mat');
X = data.Y; % SNPs
Y = data.X; % expression rate
n = data.n; % number of samples

n_iter = 1;
prediction_errors = zeros(10, 1);

% split training / test set 
training_ratio = 0.8;


for iter=1:n_iter
n_tr = floor(training_ratio * n);
n_te = n - n_tr;

rand_perm = randperm(n);
X_tr = X(rand_perm(1:n_tr), :);
Y_tr = Y(rand_perm(1:n_tr), :);
X_te = X(rand_perm(n_tr:end), :);
Y_te = Y(rand_perm(n_tr:end), :);


fprintf('sCGGM running ToyData with %d-fold cross-validation...\n', kcv);

option.verbose = true; 
option.maxiter = 500; 

opt = scggm_cv(X_tr, Y_tr, kcv, lambda1_seq, lambda2_seq, option);

% compute prediction errors
[~, e] = scggm_predict(opt.Theta, opt.intercept, X_te, Y_te);
fprintf('sCGGM demo completed, test set prediction error: %g\n', e); 
prediction_errors(iter) = e;

end
prediction_errors

