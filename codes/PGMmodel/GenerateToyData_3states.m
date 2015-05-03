clear all
close all


%% Graph parameters
p = 5; % number of cts variables
q = 7; % number of discrete variables
n_tr = 1000; % # of trainig samples
n_te = 0; % # of test samples
L = 3 * ones(q,1); % levels in each categorical variable

GenerateToyData_fn_3states(p, q, L, n_tr, n_te);

% save ToyData_v2 ToyData