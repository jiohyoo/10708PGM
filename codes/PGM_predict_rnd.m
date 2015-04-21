function [prederr] = PGM_predict_rnd(theta, alpha1, beta, betad, Y, D)

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

mu = inv(B44) * (gamma');
sigma = inv(B44);

for i = 1: size(mu,2)
    r(:,i) = mvnrnd(mu(:,i), sigma);
end

Y_hat_all = r';
res = Y - Y_hat_all;
prederr = sum(sum(res.^2)) / size(Y,2) / n;