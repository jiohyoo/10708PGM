%--------------------------------------------------------------------------
% compute covariance decomposition for gene-expression given an sCGGM, 
% SNP data and gene-expression data. 
%--------------------------------------------------------------------------

function [Cov] = scggm_cov_decompose(Theta, x, y, centered_input)

if nargin < 4
	centered_input = false;
end

iThetayy = inv(Theta.yy); 

N = size(x, 1);
if size(y, 1) ~= N
	fprintf('sCGGM:error! Genotype and expression data sample size inconsistent!\n');
	Cov = {}; 
	return;  
end

if ~centered_input
	y = y - repmat(mean(y), N, 1);
	x = x - repmat(mean(x), N, 1);
end

Cov.Overall = y'*y; 
Cov.Network_Induced = iThetayy * N; 
Cov.SNP_Induced = iThetayy * Theta.xy'* (x'*x) * Theta.xy * iThetayy; 

end
