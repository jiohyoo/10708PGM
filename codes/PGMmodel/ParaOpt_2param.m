function opt = ParaOpt_2param(ToyData, lambda1_seq, lambda2_seq, kcv, alg)


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
x = paramToVecv5_PGM(beta, betad, theta, alpha1, L, p);



%% Data split for 5-fold CV
cv_indices = crossvalind('Kfold', n_tr, kcv);
cverr = zeros(length(lambda1_seq) * length(lambda2_seq));
out_all = zeros(length(lambda1_seq) * length(lambda2_seq));
cverr2 = zeros(length(lambda1_seq), length(lambda2_seq));
out_all2 = zeros(length(lambda1_seq), length(lambda2_seq));
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

[tmp_X, tmp_Y] = meshgrid(1:length(lambda1_seq), 1:length(lambda2_seq));
index_matrix = [tmp_X(:) tmp_Y(:)];

for index_matrix_row = 1 : size(index_matrix, 1)
		lam1_idx = index_matrix(index_matrix_row, 1);
		lam2_idx = index_matrix(index_matrix_row, 2);
	    lam1 = lambda1_seq(lam1_idx);
		lam2 = lambda2_seq(lam2_idx);
		out_sum = 0;
		cverr_sum = 0;
    
    	for cv_idx = 1: kcv
        
        	X_tr_CV = X_tr(cv_indices ~= cv_idx, :);
        	D_tr_CV = D_tr(cv_indices ~= cv_idx, :);
        	X_te_CV  = X_tr(cv_indices == cv_idx, :);
        	Y_te_CV  = Y_tr(cv_indices == cv_idx, :);
        	D_te_CV  = D_tr(cv_indices == cv_idx, :);
        
	        n_tr_CV = size(X_tr_CV,1);
        
			pnopt_opts = pnopt_optimset('max_iter', 200);
        
    	    if strcmp(alg, 'PNOPT')
        	    smoothF = @(x)lhoodTfocsv5_PGM(x, D_tr_CV, X_tr_CV, L, n_tr_CV, p, q);
        	    nonsmoothH = @(varargin)tfocsProxGroupv6_PGM_2param(lam1, lam2, L, n_tr_CV, p, q, varargin{:} ); % only returns value of nonsmooth
        	    [ xopt, out, opts ] = pnopt(smoothF, nonsmoothH, x, pnopt_opts);
        	    [beta, betad, theta, alpha1] = vecToParamv5_PGM(xopt, L, n_tr_CV, p, q);
        	    
        	else
        	    %smoothF = @(x)lhoodTfocsv5_PGM(x, D_tr_CV, X_tr_CV, L, n_tr_CV, p, q);
        	    %nonsmoothH = @(varargin)tfocsProxGroupv6_PGM(lam, L, n_tr_CV, p, q, varargin{:} ); % only returns value of nonsmooth
        	    %[xopt, out, opts]=tfocs(smoothF, {}, nonsmoothH, x,opts);
        	    %[beta, betad, theta, alpha1] = vecToParamv5_PGM(xopt, L, n_tr_CV, p, q);
            
        	end
        
        	
			err = PGM_predict(theta, alpha1, beta, betad, X_te_CV, D_te_CV);
			out_sum = out_sum + out;
			cverr_sum = cverr_sum + err;

	        %disp(['CV #' num2str(cv_idx) ' - err : ' num2str(cverr(lam1_idx, lam2_idx, cv_idx))]);
		end
		cverr(index_matrix_row) = cverr_sum / kcv;
		out_all(index_matrix_row) = out_sum / kcv;

    
	    %if cverr2(lam1_idx, lam2_idx) < minerr
	    %    minerr = cverr2(lam1_idx, lam2_idx);
	    %    opt_lambda1 = lam1;
		%	opt_lambda2 = lam2
	    %end
    %end
end

for index_matrix_row = 1 : size(index_matrix, 1)
	lam1_idx = index_matrix(index_matrix_row, 1);
	lam2_idx = index_matrix(index_matrix_row, 2);
	cverr2(lam1_idx, lam2_idx) = cverr(index_matrix_row)
	out_all2(lam1_idx, lam2_idx) = out_all(index_matrix_row)
end

[opt_val1, opt_idx1] = min(cverr2);
[opt_val2, opt_idx2] = min(opt_val1);
opt_lambda1 = lambda1_seq(opt_idx1(opt_idx2));
opt_lambda2 = lambda2_seq(opt_idx2);

%% out
opt.out = out_all2;
opt.cverr = cverr;
opt.cverr2 = cverr2;
opt.minerr = minerr;
opt.lambda1 = opt_lambda1;
opt.lambda2 = opt_lambda2;

