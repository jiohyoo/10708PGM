clear all
close all


%% Graph parameters
p = 10; % number of cts variables
q = 10; % number of discrete variables
n_tr = 2000; % # of trainig samples
n_te = 200; % # of test samples
L = 2 * ones(q,1); % levels in each categorical variable

GenerateToyData_fn(p, q, L, n_tr, n_te);

% save ToyData_v2 ToyData