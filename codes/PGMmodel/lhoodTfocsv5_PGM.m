function [f, g]=lhoodTfocsv5_PGM(x,D,X,L,n,p,q )
[beta betad theta alpha1]= vecToParamv5_PGM(x,L,n,p,q);
f=flhoodv5_PGM(beta,betad,theta,alpha1,D,X,L,n,p,q);
if nargin>1
    [gradbeta gradbetad gradtheta gradalpha1]= fgradlhoodv5_PGM(beta,betad,theta,alpha1,D,X,L,n,p,q);
    g=paramToVecv5_PGM(gradbeta,gradbetad,gradtheta,gradalpha1,L,n,p,q);
end

