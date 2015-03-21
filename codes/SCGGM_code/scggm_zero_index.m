%--------------------------------------------------------------------------
% return all the nonzero values of an sCGGM estimator
% entries with absolute values <=eps will be recognized as 0
%--------------------------------------------------------------------------

function nz_theta = scggm_zero_index(theta, eps)
if nargin < 2
	eps = 0; 
end

idsxy = find(abs(theta.xy)<=eps);
idsyy = find(abs(theta.yy)<=eps);

nz_theta.xy = idsxy;
nz_theta.yy = idsyy;
