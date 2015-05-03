clear all
close all


%% Graph parameters
p = 100; % number of cts variables
q = 500; % number of discrete variables
n_tr = 50; % # of trainig samples
n_te = 20; % # of test samples
L = 2 * ones(q,1); % levels in each categorical variable

[ToyData] = GenerateToyData_fn(p, q, L, n_tr, n_te);

% save ToyData_v2 ToyData