%--------------------------------------------------------------------------
% randomly initializes sCGGM parameters Theta_xy and Theta_yy
%--------------------------------------------------------------------------

function Theta0 = scggm_initialize(J, K )

	if ( K <= 100 )
		Ai = sprandsym(K,0.01);
		Theta0.xy = sprand( J, K, 0.01 );
	else
		Ai = sprandsym(K,0.001);
		Theta0.xy = sprand( J, K, 0.001 );
	end
	Theta0.yy = 0.01*Ai*Ai' + 0.7*speye(K,K); 
end
