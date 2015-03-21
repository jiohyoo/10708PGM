%--------------------------------------------------------------------------
% calculates the penalty of sCGGM for Theta_xy and Theta_yy
% penalty = lambda1 * |Theta_xy|_1 + lambda2 * |Theta_yy|_1 
% where ||_1 is l1-norm for Theta_xy and 
% l1-norm minus the sum of diagonal entries for Theta_yy. 
%--------------------------------------------------------------------------

function p = scggm_penalty(x, lambda1, lambda2)

p = lambda1 * sum(abs(x.xy(:))) + lambda2 * sum( abs( x.yy(:) ) ) - lambda2 * sum(abs(diag(x.yy)));
