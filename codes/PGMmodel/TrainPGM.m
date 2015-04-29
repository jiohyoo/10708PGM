function opt = TrainPGM(ToyData, alg, opt_lambda)

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


opts.alg = alg;


%% for optimal lambda, train over whole train data set

if strcmp(alg, 'PNOPT')
    smoothF = @(x)lhoodTfocsv5_PGM(x, D_tr, X_tr, L, n_tr, p, q);
    nonsmoothH = @(varargin)tfocsProxGroupv6_PGM(opt_lambda, L, n_tr, p, q, varargin{:} ); % only returns value of nonsmooth
    [ xopt, out, opts ] = pnopt(smoothF, nonsmoothH, x);
    [beta, betad, theta, alpha1] = vecToParamv5_PGM(xopt, L, n_tr, p, q);
    
else
    smoothF = @(x)lhoodTfocsv5_PGM(x, D_tr, X_tr, L, n_tr, p, q);
    nonsmoothH = @(varargin)tfocsProxGroupv6_PGM(opt_lambda, L, n_tr, p, q, varargin{:} ); % only returns value of nonsmooth
    [xopt, out, opts]=tfocs(smoothF, {}, nonsmoothH, x,opts);
    [beta, betad, theta, alpha1] = vecToParamv5_PGM(xopt, L, n_tr, p, q);
    
end


%% out
opt.beta = beta;
opt.betad = betad;
opt.theta = theta;
opt.alpha1 = alpha1;
opt.lambda = opt_lambda;
opt.out = out;
