%--------------------------------------------------------------------------
% perform inference to obtain indirect SNP perturbations given an sCGGM
%--------------------------------------------------------------------------

function [Beta] = scggm_indirect_SNP_overall(Theta)
Beta = - Theta.xy * inv(Theta.yy); 


end
