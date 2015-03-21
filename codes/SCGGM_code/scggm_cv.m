%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Estimates a sparse CGGM with cross validation to select optimal
%  regularization parameters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ OPT ] = scggm_cv( x, y, kcv, lambda1_seq, lambda2_seq, option)
default_lambdaseq = [0.32, 0.16, 0.08, 0.04, 0.02, 0.01]; 
if ( nargin <4 || nargin == 5)
    maxiter = 1000;
    tol     = 1e-7;
    verbose = false;
    eta     = 1.5; 
    centered_input = false; 
    ifrefit = true; 
    Theta0  = scggm_initialize(size(x, 2), size(y, 2));
    if nargin < 4
	lambda1_seq = default_lambdaseq; 
	lambda2_seq = default_lambdaseq; 
    elseif ~isa(lambda1_seq, 'double') || ~isa(lambda2_seq, 'double')
	fprintf('sCGGM: error!lambda1_seq lambda2_seq must be vectors\n');
	OPT = {}; return;
    end
else
    if nargin == 4
	if isa(lambda1_seq, 'struct')
		option = lambda1_seq; 
		lambda1_seq = default_lambdaseq; 
		lambda2_seq = default_lambdaseq; 
	else
		fprintf('sCGGM: error! scggm_cv must accept lambda1_seq lambda2_seq simultaneously!\n');
		OPT = {}; return; 
	end
    end
    if  isfield(option,'maxiter')
        maxiter = option.maxiter; 
    else
        maxiter = 1000;
    end
    if isfield(option, 'centered_input')
        centered_input = option.centered_input; 
    else
	centered_input = false; 
    end
    if isfield(option, 'ifrefit')
	ifrefit = option.ifrefit; 
    else
	ifrefit = true; 
    end
    if isfield(option, 'Theta0')
        Theta0 = option.Theta0; 
    else
        Theta0 = scggm_initialize(size(x, 2), size(y, 2));
    end
    if isfield(option, 'tol')
        tol = option.tol; 
    else
        tol = 1e-7;
    end
    if isfield(option, 'verbose')
        verbose = option.verbose; 
    else
        verbose = false;
    end
    if isfield(option, 'eta')
        eta  	= option.eta; 
    else
        eta     = 1.5; 
    end
end

%% center the input data if not centered 
N0 = size(x, 1);
if size(y, 1) ~= N0
	fprintf('sCGGM:error! Input data sample size inconsistent!\n');
	OPT = {}; return;  
end
if kcv < 3 || kcv >= N0
	fprintf('sCGGM:error! Cross validation cannot be %d fold \n',kcv); 
	OPT = {}; return; 
end
if ~centered_input
	y0 = y - repmat(mean(y), N0, 1);
	x0 = x - repmat(mean(x), N0, 1);
else
	x0 = x; 
	y0 = y; 
end
J  = size(x, 2);
K  = size(y, 2);


% cross validation index
cv_indices  = crossvalind('Kfold', N0, kcv);
cverr       = zeros(length(lambda1_seq), length(lambda2_seq));
minerr      = 1e99;

if verbose 
	fprintf('J = %d, K = %d, sample size = %d\n',J, K, N0); 
end

%% k-fold cross validation 
for i = 1:length(lambda1_seq)
    for j = 1:length(lambda2_seq)
        lambda1 = lambda1_seq(i); 
        lambda2 = lambda2_seq(j);
        
        for ff = 1:kcv
		% extract centered and uncentered training data
		cxtr = x0(cv_indices ~= ff, :); 
		cytr = y0(cv_indices ~= ff, :); 
		xtr  = x(cv_indices ~= ff, :); 
		ytr  = y(cv_indices ~= ff, :); 

		% extract uncentered cross-validation data
		xcv  = x(cv_indices == ff, :);
		ycv  = y(cv_indices == ff, :); 
		
		% estimate a sparse estimate of Theta_xy and Theta_yy
		raw_Theta = scggm_sparse_step(lambda1, lambda2,cxtr, cytr, maxiter, tol, verbose, eta, Theta0);
 
		if ifrefit
			% refit the parameters
			zero_theta = scggm_zero_index(raw_Theta);
			Theta = scggm_refit_step(cxtr, cytr, zero_theta, maxiter, tol, verbose, eta, raw_Theta); 
            	else
			Theta = raw_Theta; 
		end

		% compute cross-validation error
		intercept = mean(ytr) + mean(xtr) * (Theta.xy * inv(Theta.yy));
            	[Eycv, e] = scggm_predict(Theta, intercept, xcv, ycv); 
		cverr(i,j) = cverr(i,j) + e; 
        end
        
        cverr(i,j) = cverr(i,j) / kcv; 
        if cverr(i,j) < minerr
            minerr      = cverr(i,j); 
            opt_lambda1 = lambda1; 
            opt_lambda2 = lambda2; 
        end
	if verbose     	   
            fprintf('lambda_1 = %.2f\t lambda_2 = %.2f\t cross validation error = %.3f\n', lambda1, lambda2, cverr(i,j)); 
	end
    end
end

if verbose 
    fprintf('\ntraining sCGGM with optimal regularization parameters: \n'); 
    fprintf('optimal lambda_1 = %.2f\t optimal lambda_2 = %.2f... \n', opt_lambda1, opt_lambda2); 
end

[raw_Theta] = scggm_sparse_step( opt_lambda1, opt_lambda2, x0, y0,  maxiter, tol, verbose, eta, Theta0);

if ifrefit
	zero_theta = scggm_zero_index(raw_Theta);
	Theta = scggm_refit_step(x0, y0, zero_theta, maxiter, tol, verbose, eta, raw_Theta);
else
	Theta = raw_Theta; 
end

OPT.Theta 	= Theta; 
OPT.lambdas     = [ opt_lambda1, opt_lambda2 ];
OPT.intercept	= mean(y) + mean(x) * (Theta.xy * inv(Theta.yy)); 

end
