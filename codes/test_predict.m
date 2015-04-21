function [prederr] = test_predict(theta, alpha1, beta, betad, Y, D)

% n = size(Y,1);
% Dtheta = D * theta;
% gamma = (Dtheta + repmat(alpha1',n,1));
% B44 = -beta + diag(betad);
% % Y_hat_all = inv(B44) * (gamma');
% Y_hat_all = (B44) \ (gamma');
% res = Y - Y_hat_all';
% prederr = sum(sum(res.^2)) / size(Y,2) / n;

n = size(Y,1);
Dtheta = D * theta;
gamma = (Dtheta + repmat(alpha1',n,1));
B44 = -beta + diag(betad);
for i = 1: size(B44,1)
    for j = 1: size(B44,2)
        if (i>j)
            B44(i,j) = B44(j,i);
        end
    end
end

% Y_hat_all = inv(B44) * (gamma');
mu = (B44) \ (gamma');
sigma = sqrt(inv(B44));

for i = 1: size(mu,2)
    for j = 1: n_sample
        R( = mvnrnd(mu(:,i),sigma);
    end
end

res = Y - Y_hat_all';
prederr = sum(sum(res.^2)) / size(Y,2) / n;