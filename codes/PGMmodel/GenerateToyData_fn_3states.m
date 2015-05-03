function GenerateToyData_fn_3states(p, q, L, n_tr, n_te)

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
gen_done = false;
while gen_done == false
    thcts = rand(p, p) * 2 - 1;
    thcts = thcts * thcts';
    cts_edge_prob = 0.6;

    for i = 1:p
        for j = i+1:p
            if rand() > cts_edge_prob
                thcts(i, j) = 0;
                thcts(j, i) = 0;
            end
        end
    end

    if sum(eig(thcts) <= 0) == 0
        gen_done = true;
    end
end


thdiscts = cell(p,q); % theta or rho in the paper
maskDisCts = sparse(p,q); % ground-truth relation indicator
discts_edge_prob = 0.4;

for i = 1:p
    for j=1:q
        if rand() < discts_edge_prob
            maskDisCts(i,j) = 1;
            if rand() < 0.5
                %%%%%%%%%%%%%%%%%%%%% for each of three states
                thdiscts{i,j} = [1.0, 1.2, 3.5];
            else
                thdiscts{i,j} = [-3.5, -3.1, 1.0];
            end
        end
    end
end

thdis = cell(q,q); % phi in the paper - we don't care about this term
maskDis = sparse(diag(ones(q-1,1),1));
%maskDis = sparse(diag(ones(q,1),1));
[R J] = find(maskDis);
for e = 1 : length(R)
        thdis{R(e),J(e)}= .5 * [1 0 -1 ; 0 1 0 ; -1 0 1]; % attractive %%%%%%%%%%%%%%%%%%%%%
end
%% sample
mixSample;

save ToyData4states X Y D p q L n thcts thdiscts maskDisCts maskDis
% ToyData.Y_tr = Y;
% ToyData.X_tr = X;
% ToyData.D_tr = D;
% ToyData.p = p;
% ToyData.q = q;
% ToyData.L = L;
% ToyData.n_tr = n;
% ToyData.thcts = thcts;
% ToyData.maskDisCts = maskDisCts;
% ToyData.maskDis = maskDis;
% ToyData.nodePot = nodePot;
% ToyData.edgePot = edgePot;
% ToyData.edgeStruct = edgeStruct;
% ToyData.burnIn = burnIn;
% ToyData.thdiscts = thdiscts;
% 
% % save ToyData_v2 ToyData
% 
% %% generate additional data for test set
% % clear all
% % load ToyData_v2
% 
% % ToyData.n_te = 200;
% ToyData.n_te = n_te;
% 
% [X_te, Y_te, D_te] = GenerateTestSample(ToyData);
% ToyData.X_te = X_te;
% ToyData.Y_te = Y_te;
% ToyData.D_te = D_te;
% 
% % save ToyData_v2 ToyData
% 
