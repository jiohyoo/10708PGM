function [f, g]=lhoodTfocsv5_PGM(x,D,X,Y,L,n,p,q )
[beta betad theta phi alpha1 alpha2]= vecToParamv5(x,L,n,p,q);
f=flhoodv5_PGM(beta,betad,theta,phi,alpha1,alpha2,D,X,Y,L,n,p,q);
if nargin>1
    [gradbeta gradbetad gradtheta gradphi gradalpha1 gradalpha2  ]= fgradlhoodv5_PGM(beta,betad,theta,phi,alpha1,alpha2,D,X,Y,L,n,p,q);
    g=paramToVecv5(gradbeta,gradbetad,gradtheta,gradphi,gradalpha1,gradalpha2,L,n,p,q);
end

