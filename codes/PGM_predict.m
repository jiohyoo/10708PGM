function [prederr] = PGM_predict(theta, alpha1, beta, betad, Y, D)

n = size(Y,1);
Dtheta = D * theta;
gamma = (Dtheta + repmat(alpha1',n,1));
B44 = -beta + diag(betad);

for i = 1: 10
    for j = 1: 10
        if i>j
           B44(i,j) = B44(j,i); 
        end
    end
end

% Y_hat_all = inv(B44) * (gamma');
Y_hat_all = (B44) \ (gamma');
res = Y - Y_hat_all';
prederr = sum(sum(res.^2)) / size(Y,2) / n;