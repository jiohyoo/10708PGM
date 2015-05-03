clear all
close all


%% Graph parameters
p = 10; % number of cts variables
q = 10%60; % number of discrete variables
n_tr = 1000;%1000%15; % # of trainig samples
n_te = 2; % # of test samples
L = 2 * ones(q,1); % levels in each categorical variable


load maskDisCts_ex maskDisCts
GenerateToyData_fn3(p, q, L, n_tr, n_te, maskDisCts);

% save ToyData_v2 ToyData