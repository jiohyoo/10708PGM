function opt = TrainPGM_jason(ToyData, alg, opt_lambda)

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

n = n_tr;
D = D_tr;
X = X_tr;
Y = Y_tr;

%% Opt alg
opts.alg = alg;


%% Init opt variables
theta = zeros(Ltot,p); % cts-dis params
beta = zeros(p,p); % negative of the precision matrix
betad = ones(p,1); % diagonal of the precision matrix
alpha1 = zeros(p,1); % cts node potential param
alpha2=zeros(Ltot,1); % dis node potential param
phi=zeros(Ltot,Ltot); % dis edge potential params
Lsum = [0; cumsum(L)];
x=paramToVecv5(beta,betad,theta,phi,alpha1,alpha2,L,n,p,q);


options.tfocs_opts.max_iter = 10;

%% for optimal lambda, train over whole train data set


smoothF= @(x)lhoodTfocsv5(x,D,X,Y,L,n,p,q);
nonsmoothH=@(varargin) tfocsProxGroupv6(opt_lambda,L,n,p,q, varargin{:} ); % only returns value of nonsmooth
opts.alg='N83';  opts.maxIts=800; opts.printEvery=100; opts.saveHist=true;
opts.restart=-10^4;
opts.tol=1e-10;
% [xopt out opts]=tfocs(smoothF, {}, nonsmoothH, x,opts);
[ xopt, out, opts ] = pnopt( smoothF, nonsmoothH, x);
[beta betad theta phi alpha1 alpha2]= vecToParamv5(xopt,L,n,p,q);



%% out
opt.beta = beta;
opt.betad = betad;
opt.theta = theta;
opt.alpha1 = alpha1;
opt.lambda = opt_lambda;
opt.out = out;




