function opt = syntheticExp2_PGM_CV_fn(ToyData, lambda_seq, kcv)


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
% alpha2=zeros(Ltot,1); % dis node potential param
% phi=zeros(Ltot,Ltot); % dis edge potential params
Lsum = [0; cumsum(L)];
x = paramToVecv5_PGM(beta, betad, theta, alpha1, L, n_tr, p, q);



%% Data split for 5-fold CV

cv_indices = crossvalind('Kfold', n_tr, kcv);
cverr = zeros(1, length(lambda_seq));
minerr = 1e99;

for lam_idx = 1 : length(lambda_seq)
    lam = lambda_seq(lam_idx);
    
    for cv_idx = 1: kcv
         
		X_tr_CV = X_tr(cv_indices ~= cv_idx, :); 
		D_tr_CV = D_tr(cv_indices ~= cv_idx, :); 
		X_te_CV  = X_tr(cv_indices == cv_idx, :);
        Y_te_CV  = Y_tr(cv_indices == cv_idx, :); 
		D_te_CV  = D_tr(cv_indices == cv_idx, :); 
        
        n_tr_CV = size(X_tr_CV,1);
        
        smoothF = @(x)lhoodTfocsv5_PGM(x, D_tr_CV, X_tr_CV, L, n_tr_CV, p, q);
        nonsmoothH = @(varargin)tfocsProxGroupv6_PGM(lam, L, n_tr_CV, p, q, varargin{:} ); % only returns value of nonsmooth
        % opts.alg='N83';  opts.maxIts=800; opts.printEvery=100; opts.saveHist=true;
        % opts.restart=-10^4;
        % opts.tol=1e-10;
        % [xopt out opts]=tfocs(smoothF, {}, nonsmoothH, x,opts);
        % [ xopt, out, opts ] = pnopt( smoothF, nonsmoothH, x, opts );
        [ xopt, out, opts ] = pnopt(smoothF, nonsmoothH, x);
        [beta, betad, theta, alpha1] = vecToParamv5_PGM(xopt, L, n_tr_CV, p, q);
        
        
        
        %% eval lambda
%         [prederr_tmp] = PGM_predict(theta, alpha1, beta, betad, Y_te_CV, D_te_CV);
            prederr_tmp = out;
%         fprintf('PGM_model - lambda: %g\n, prediction error: %g\n', lam, prederr_tmp);
        
        cverr(lam_idx) = cverr(lam_idx) + prederr_tmp; 
        
    end
    
    
    cverr(lam_idx) = cverr(lam_idx) / kcv;
    if cverr(lam_idx) < minerr
        minerr = cverr(lam_idx);
        opt_lambda = lam;
    end
    
end







%% for optimal lambda, 

smoothF = @(x)lhoodTfocsv5_PGM(x, D_tr, X_tr, L, n_tr, p, q);
nonsmoothH = @(varargin)tfocsProxGroupv6_PGM(opt_lambda, L, n_tr, p, q, varargin{:} ); % only returns value of nonsmooth

[xopt, out, opts] = pnopt(smoothF, nonsmoothH, x);
[beta, betad, theta, alpha1] = vecToParamv5_PGM(xopt, L, n_tr, p, q);





%% out
opt.beta = beta;
opt.betad = betad;
opt.theta = theta;
opt.alpha1 = alpha1;
opt.lambda = opt_lambda;



% [prederr_fin, opt_lambda_idx] = min(prederr_all);
% opt_lambda = lambda_seq(opt_lambda_idx);
% fprintf('---------------------------------------\n');
% fprintf('PGM_model - opt lambda: %g\n, fin prediction error: %g\n', opt_lambda, prederr_fin);