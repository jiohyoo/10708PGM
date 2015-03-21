%--------------------------------------------------------------------------
% evaluate sCGGM objective and/or gradient
%--------------------------------------------------------------------------

function [ value, flag, grad] = scggm_evaluate(theta, Sx, Sxy, Sy, N, gradient, verbose)

flag        = 0;
[ cyy, p ]  = chol(theta.yy);

if ( p > 0 )
	if strcmp(gradient, 'y') == 1 && verbose
		fprintf( 'sCGGM: Theta_yy not positive definite!\n' );
	end
	flag        = 1;
	value       = inf;
	grad        = theta;
	return;
end

logdetyy = 2 * sum(log(diag(cyy) ));

if ( isnan(logdetyy) || isinf(logdetyy) )
	if verbose
		fprintf( 'sCGGM: logdet Theta_yy is Nan or Inf!\n' );
	end
	flag = 1;
	value = inf;
	grad = theta;
	return; 
end

icyy	 = cyy \ eye(size(cyy,2));
ithetayy = icyy * icyy';
txyityy  = theta.xy*ithetayy;
XtXth    = Sx*txyityy;
txyXtXth = theta.xy'*Sx*txyityy;

l1 = trace( theta.yy*Sy );
l2 = trace( Sxy*theta.xy' );
l3 = trace( txyXtXth );
value = 0.5*l1 + l2 + 0.5*l3 - 0.5*N*logdetyy ;
value = value / N;

if strcmp('y',gradient) ==1
	grad.xy = (Sxy + XtXth)/N;
	grad.yy = 0.5*(Sy - N*ithetayy - ithetayy*txyXtXth)/N;
else
	grad = [];
end

