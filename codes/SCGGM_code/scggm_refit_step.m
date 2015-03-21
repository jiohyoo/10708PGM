%--------------------------------------------------------------------------
% refit sCGGM estimators returned by function scggm_sparse_step()
%
% the parameter is estimated with the constraint that
% at positions z_theta.xy and z_theta.yy, parameters are 0
% use for refitting after scggm_sparse_step() procedure
%--------------------------------------------------------------------------

function [Theta, obj] = scggm_refit_step( cx, cy, z_theta, maxiter, tol, verbose, eta, Theta0)
Sx	= cx'*cx; 
Sy 	= cy'*cy; 
Sxy 	= cx'*cy; 
N 	= size(cx, 1);

theta = Theta0;
bconv = 0;
theta.xy(z_theta.xy) = 0; 
theta.yy(z_theta.yy) = 0; % constrain the sparsity pattern of the variable
obj  = zeros(maxiter, 1);
L    = 1;
nobj = 10;

thk_0 	= 2/3; 
ls_maxiter = 300; 

[obj1, init_flag]  = scggm_evaluate( theta, Sx, Sxy, Sy, N, 'n', verbose);
if init_flag == 1 && verbose 
	fprintf('sCGGM: error! refit initial Theta_yy not positive definite!\n');
end

obj(1)  = obj1; 
xk      = theta;
zk      = theta;
thk     = thk_0;

for iter = 2:maxiter
	% thk  = 2*thk / (2+ thk);
	thk  = (sqrt( thk^4 + 4 * thk^2 ) - thk^2) / 2; % momentum of acceleration
	y.xy = (1-thk)*xk.xy + thk*zk.xy;
	y.yy = (1-thk)*xk.yy + thk*zk.yy;
	[ fyk, flagy, grady] = scggm_evaluate( y, Sx, Sxy, Sy, N, 'y', verbose); % compute the objective and gradient for y
	grady.xy(z_theta.xy) = 0;
	grady.yy(z_theta.yy) = 0;

 	% line search
	ik = 0;
	while ( 1 )
        	% gradient descent
		zk_grady.xy = y.xy - 1/(thk*L) * grady.xy;
		zk_grady.yy = y.yy - 1/(thk*L) * grady.yy;
		zk1 = zk_grady;

        	% gradient descent
		xk1.xy = y.xy - 1/L *grady.xy;
		xk1.yy = y.yy - 1/L *grady.yy;
		[fxk1, flagxk1] = scggm_evaluate(xk1, Sx, Sxy, Sy, N , 'n', verbose);      
		[~, flagzk1] = chol(zk1.yy);
                
		if ( flagzk1 == 0 && flagy ==0 && flagxk1 ==0 ) % xk1,zk1,y all positive definite
			xk1_y.xy = xk1.xy - y.xy;
			xk1_y.yy = xk1.yy - y.yy;
			lfxk1_y  = fyk + grady.xy(:)' * (xk1_y.xy(:)) + grady.yy(:)'*(xk1_y.yy(:));
			diffxk1.xy = xk1.xy - y.xy;
			diffxk1.yy = xk1.yy - y.yy;
			RHS = lfxk1_y +  L/2 * (sum(diffxk1.xy(:).^2) + sum(diffxk1.yy(:).^2));
			LHS = fxk1;
			if ( LHS <= RHS + tol)
				xk = xk1;
				zk = zk1;
				bconv = 1;
				break; % line search converged
			end
		end
                
		ik = ik + 1;
		if ( ik > ls_maxiter )
			if verbose
				fprintf( 'sCGGM: refit line search not converging,ik = %d\n',ik);
			end
			bconv = 0;
			iter  = max(1, iter - 1); 
			Theta = xk;
			break;
		end
                
        	L = L*eta;
	end
	obj(iter) = fxk1;
	if bconv ==0
		break; 
	end
	if ( iter > nobj + 1 )
		value           = obj( iter );
		prevVals        = obj( iter - nobj);
		avgimprovement  = abs( prevVals - value )/ nobj;
		relAvgImpr      = avgimprovement / abs( value );
		if ( relAvgImpr < tol )
			bconv = 1;
			break;
		end
	end
end

Theta = xk;
obj   = obj(1:iter);
