%--------------------------------------------------------------------------
% soft threshold operation on Theta_xy and Theta_yy
% thresholded at value c1 (for Theta_xy) and c2 (for Theta_yy)
%--------------------------------------------------------------------------
function B = scggm_soft_threshold(theta, c1, c2)

    Bxy=zeros(size(theta.xy));
    pos_idx=theta.xy(:)>c1;
    neg_idx=theta.xy(:)<-c1;
    
    Bxy(pos_idx)=theta.xy(pos_idx)-c1;
    Bxy(neg_idx)=theta.xy(neg_idx)+c1;

    Byy = diag(diag(theta.yy));
    uyy = triu( theta.yy,1 );
    pos_idx=uyy(:)>c2;
    neg_idx=uyy(:)<-c2;

    Byy(pos_idx) = uyy(pos_idx)-c2;
    Byy(neg_idx) = uyy(neg_idx)+c2;
    Byy = Byy + Byy'-diag(diag(theta.yy));

    B.xy = (Bxy);
    B.yy = (Byy);

end
