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

%% Opt alg
opts.alg = alg;


%% Init opt variables
theta = zeros(Ltot,p); % cts-dis params
beta = zeros(p,p); % negative of the precision matrix
betad = ones(p,1); % diagonal of the precision matrix
alpha1 = zeros(p,1); % cts node potential param
% alpha2=zeros(Ltot,1); % dis node potential param
% phi=zeros(Ltot,Ltot); % dis edge potential params
Lsum = [0; cumsum(L)];
x = paramToVecv5_PGM(beta, betad, theta, alpha1, L, p);

options.tfocs_opts.max_iter = 10;

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
