function opt = ParaOpt_jason(ToyData, lambda_seq, kcv, alg)


%% Data
X_tr = ToyData.X_tr;
Y_tr = ToyData.Y_tr;
D_tr = ToyData.D_tr;
X_te = ToyData.X_te;
Y_te = ToyData.Y_te;
D_te = ToyData.D_te;
p = ToyData.p;
q = ToyData.q;
L = ToyData.L;
n_tr = ToyData.n_tr;

e = ones(n_tr,1);
Ltot = sum(L);



%% Init opt variables
theta = zeros(Ltot,p); % cts-dis params
beta = zeros(p,p); % negative of the precision matrix
betad = ones(p,1); % diagonal of the precision matrix
alpha1 = zeros(p,1); % cts node potential param
alpha2=zeros(Ltot,1); % dis node potential param
phi=zeros(Ltot,Ltot); % dis edge potential params
Lsum = [0; cumsum(L)];
x=paramToVecv5(beta,betad,theta,phi,alpha1,alpha2,L,n_tr,p,q);



%% Data split for 5-fold CV
cv_indices = crossvalind('Kfold', n_tr, kcv);
cverr = zeros(length(lambda_seq), kcv);
cverr2 = zeros(length(lambda_seq), 1);
minerr = 1e99;
testerr = [];


%% define settings for optimizaiton algs
% opts.maxIts = 800;
% opts.printEvery = 100;
% opts.saveHist = true;
% opts.restart = -10^4;
% opts.tol = 1e-10;
% opt_lambda = zeros(1, n_opt_algs);
opts.alg = alg;

%%


for lam_idx = 1 : length(lambda_seq)
    lam = lambda_seq(lam_idx);
    
    for cv_idx = 1: kcv
        
        X_tr_CV = X_tr(cv_indices ~= cv_idx, :);
        Y_tr_CV  = Y_tr(cv_indices ~= cv_idx, :);
        D_tr_CV = D_tr(cv_indices ~= cv_idx, :);
        X_te_CV  = X_tr(cv_indices == cv_idx, :);
        Y_te_CV  = Y_tr(cv_indices == cv_idx, :);
        D_te_CV  = D_tr(cv_indices == cv_idx, :);
        
        n_tr_CV = size(X_tr_CV,1);
        
        %% Init opt variables
        theta = zeros(Ltot,p); % cts-dis params
        beta = zeros(p,p); % negative of the precision matrix
        betad = ones(p,1); % diagonal of the precision matrix
        alpha1 = zeros(p,1); % cts node potential param
        alpha2=zeros(Ltot,1); % dis node potential param
        phi=zeros(Ltot,Ltot); % dis edge potential params
        Lsum = [0; cumsum(L)];
        x=paramToVecv5(beta,betad,theta,phi,alpha1,alpha2,L,n_tr_CV,p,q);
        
        smoothF= @(x)lhoodTfocsv5(x,D_tr_CV,X_tr_CV,Y_tr_CV,L,n_tr_CV,p,q);
        nonsmoothH=@(varargin) tfocsProxGroupv6(lam,L,n_tr_CV,p,q, varargin{:} ); % only returns value of nonsmooth
        opts.alg='N83';  opts.maxIts=800; opts.printEvery=100; opts.saveHist=true;
        opts.restart=-10^4;
        opts.tol=1e-10;
        % [xopt out opts]=tfocs(smoothF, {}, nonsmoothH, x,opts);
        [ xopt, out, opts ] = pnopt( smoothF, nonsmoothH, x);
        [beta betad theta phi alpha1 alpha2]= vecToParamv5(xopt,L,n_tr_CV,p,q);
        
        out_all(lam_idx, cv_idx) = out;
        %% eval lambda
        cverr(lam_idx, cv_idx) = PGM_predict(theta, alpha1, beta, betad, X_te_CV, D_te_CV);
        disp(['CV #' num2str(cv_idx) ' - err : ' num2str(cverr(lam_idx, cv_idx))]);
        clear x
    end
    
    cverr2(lam_idx) = sum(cverr(lam_idx, : )) / kcv;
    
    if cverr2(lam_idx) < minerr
        minerr = cverr2(lam_idx);
        opt_lambda = lam;
    end
    
    
end



%% out
opt.out = out_all;
opt.cverr = cverr;
opt.cverr2 = cverr2;
opt.minerr = minerr;
opt.lambda = opt_lambda;

