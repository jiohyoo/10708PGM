%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The main sparse CGGM function that parses the parameter settings for the 
%  accelerated proximal gradient algorithm, estimate sCGGM parameters with
%  regularization and refit the parameters given the sparsity pattern.
%  Refit step can be skipped to focus on finding sparsity pattern of 
%  Theta_xy and Theta_yy. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OPT =  scggm( x, y, lambda_1, lambda_2, option)
%% setting algorithm running options
if ( nargin < 5)
    maxiter = 1000;
    tol     = 1e-7;
    verbose = false;
    eta     = 1.5; 
    centered_input = false; 
    ifrefit = true; 
    Theta0  = scggm_initialize(size(x, 2), size(y, 2));
else
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
    if isfield(option, 'Theta0')
	Theta0 = option.Theta0; 
    else
	Theta0 = scggm_initialize(size(x, 2), size(y, 2));
    end
end

N = size(x, 1);
if size(y, 1) ~= N
	fprintf('sCGGM:error! Input data sample size inconsistent!\n');
	OPT = {}; 
	return;  
end
%% center the input data if not centered 
if ~centered_input
	cy = y - repmat(mean(y), N, 1);
	cx = x - repmat(mean(x), N, 1);
else
	cx = x; 
	cy = y; 
end

%% estimate a sparse estimate of Theta_xy and Theta_yy
[raw_Theta] = scggm_sparse_step( lambda_1 , lambda_2,cx, cy,  maxiter, tol, verbose, eta, Theta0);
 
if ifrefit
	%% refit the parameters
	zero_theta = scggm_zero_index(raw_Theta);
	Theta = scggm_refit_step(cx, cy, zero_theta, maxiter, tol, verbose, eta, raw_Theta); 
else
	Theta = raw_Theta; 
end

%% record results 
OPT.Theta = Theta;
OPT.intercept = mean(y) + mean(x) * (Theta.xy * inv(Theta.yy)); 
