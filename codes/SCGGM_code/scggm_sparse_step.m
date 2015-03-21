%--------------------------------------------------------------------------
% estimate a sparse CGGM
%--------------------------------------------------------------------------

function [Theta, obj] = scggm_sparse_step( lambda1, lambda2, cx, cy, maxiter, tol, verbose, eta, Theta0)

Sx	= cx'*cx; 
Sy 	= cy'*cy; 
Sxy	= cx'*cy;
N 	= size(cx, 1);

nobj	= 10;
bconv	= 0;
obj	= zeros(maxiter, 1);
theta	= Theta0; 
L	= 1; 
thk_0 	= 2/3; 
ls_maxiter = 300; 

[obj1, init_flag]  = scggm_evaluate( theta, Sx, Sxy, Sy, N, 'n', verbose);
if init_flag == 1 && verbose == true
	fprintf('sCGGM: error! initial Theta_yy not positive definite!\n');
end
obj(1) = obj1 + scggm_penalty(theta, lambda1, lambda2);

xk      = theta;
zk   	= theta;
thk     = thk_0;

for iter = 2:maxiter
	thk  = (sqrt( thk^4 + 4 * thk^2 ) - thk^2) / 2; % momentum of acceleration
	y.xy = (1 - thk) * xk.xy + thk * zk.xy;
	y.yy = (1 - thk) * xk.yy + thk * zk.yy;
	[ fyk, flagy, grady] = scggm_evaluate( y, Sx, Sxy, Sy, N, 'y', verbose); % compute the objective and gradient for y
    
	% line search
	ik = 0; 
	while ( 1 )
 	       % gradient step
		zk_grady.xy = zk.xy - 1/(L*thk) * grady.xy;	
		zk_grady.yy = zk.yy - 1/(L*thk) * grady.yy;
        	% proximal step
		zk1 		= scggm_soft_threshold( zk_grady, 2*lambda1/(thk*L), 2*lambda2/(thk*L)) ;
	        % gradient step
		y_grady.xy	= y.xy - 1/L * grady.xy;
		y_grady.yy	= y.yy - 1/L * grady.yy;
       		% proximal step
        	xk1         = scggm_soft_threshold( y_grady, 2*lambda1/(L), 2*lambda2/(L));
        
		[fxk1, flagxk1] = scggm_evaluate(xk1, Sx, Sxy, Sy, N ,'n', verbose);
		[~, flagzk1]    = chol(zk1.yy);

		if ( flagzk1 == 0 && flagy ==0 && flagxk1 ==0 ) % xk1,zk1,y all positive definite
			xk1_y.xy    = xk1.xy - y.xy;
			xk1_y.yy    = xk1.yy - y.yy;	
			lfxk1_y     = fyk + grady.xy(:)'* (xk1_y.xy(:)) + grady.yy(:)'*(xk1_y.yy(:));
			diffxk1y.xy = xk1.xy - y.xy;
			diffxk1y.yy = xk1.yy - y.yy;
			RHS         = lfxk1_y + L/2 *(sum(diffxk1y.xy(:).^2) + sum(diffxk1y.yy(:).^2));
			if fxk1 <= RHS + tol
				xk = xk1;
				zk = zk1;
				bconv = 1;
				break; % line search converged
			end
        	end
        
		ik = ik + 1;
        
		if ( ik > ls_maxiter )
			if verbose
				fprintf( 'sCGGM: line search not converging,ik = %d\n',ik); 
			end
			bconv = 0;
			iter  = max(1, iter - 1); 
			Theta = xk; 
			break;
		end
		L = L * eta;
	end 
	obj(iter)  = fxk1 + scggm_penalty( xk, lambda1, lambda2);
	if bconv == 0
		break;
	end
    
	if ( iter > nobj + 1)
		value           = obj(iter);
		prevVals        = obj(iter - nobj);
		avgimprovement  = abs( prevVals - value )/nobj;
		relAvgImpr      = avgimprovement / abs( value ) ; % relative average improvement
            
		if ( relAvgImpr < tol )
			bconv = 1;
			break;
		end
	end
end  

Theta = xk; 
obj   = obj(1:iter);
