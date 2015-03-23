function [ToyData] = GenerateToyData_fn(p, q, L, n_tr, n_te)

addpath(genpath('./UGM_2009')); 
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));



%% Graph parameters
n = n_tr;
% p = 10; % number of cts variables
% q = 10; % number of discrete variables
% n = 2000; % sample size
e = ones(n,1);
% L = 2*ones(q,1); % levels in each categorical variable
Ltot = sum(L);



%% Create params
rhocts = 0.25;
thcts = eye(p);
for i = 1 : p - 1
    thcts(i, i+1) = rhocts;
    thcts(i+1, i) = rhocts;
end

thdiscts = cell(p,q);
maskDisCts = sparse(p,q);
for i = 1 : min(p,q)
    rhocd = .5;
    maskDisCts(i,i) = 1;
    thdiscts{i,i} = [rhocd; -rhocd];
end

thdis = cell(q,q);
maskDis = sparse(diag(ones(q-1,1),1));
[R J] = find(maskDis);
for e = 1 : length(R)
        thdis{R(e),J(e)}= .5 * [1 -1; -1 1]; % attractive
end
%% sample
mixSample;


ToyData.Y_tr = Y;
ToyData.X_tr = X;
ToyData.D_tr = D;
ToyData.p = p;
ToyData.q = q;
ToyData.L = L;
ToyData.n_tr = n;
ToyData.thcts = thcts;
ToyData.maskDisCts = maskDisCts;
ToyData.maskDis = maskDis;
ToyData.nodePot = nodePot;
ToyData.edgePot = edgePot;
ToyData.edgeStruct = edgeStruct;
ToyData.burnIn = burnIn;
ToyData.thdiscts = thdiscts;

% save ToyData_v2 ToyData

%% generate additional data for test set
% clear all
% load ToyData_v2

% ToyData.n_te = 200;
ToyData.n_te = n_te;

[X_te, Y_te, D_te] = GenerateTestSample(ToyData);
ToyData.X_te = X_te;
ToyData.Y_te = Y_te;
ToyData.D_te = D_te;

% save ToyData_v2 ToyData

