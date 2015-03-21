%--------------------------------------------------------------------------
% compute decomposition of the indirect SNP perturbations for a particular gene
% given an sCGGM
%--------------------------------------------------------------------------

function [Beta_k] = scggm_indirect_SNP_decompose(Theta, k)

iThetayy = inv(Theta.yy); 
Beta_k = - Theta.xy(:,k) * iThetayy(:,k)';

end
