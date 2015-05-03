
theta = opt.theta;
alpha1 = opt.alpha1;
beta = opt.beta;
betad = opt.betad;
Y = ToyData.X_te;
D = ToyData.D_te;

n = size(Y,1);
Dtheta = D * theta;
gamma = (Dtheta + repmat(alpha1',n,1));
B44 = -beta + diag(betad);

for i = 1: size(B44,1)
    for j = 1:size(B44,1)
        if i < j
            B44(j,i) = B44(i,j);
        end
    end
end

% Y_hat_all = inv(B44) * (gamma');
Y_hat_all = (B44) \ (gamma');
res = Y - Y_hat_all';
prederr = sum(sum(res.^2)) / size(Y,2) / n;

cmax = max(max(Y));
cmin = min(min(Y));
figure;imagesc(Y);colorbar;
figure;imagesc(Y_hat_all');caxis([cmin cmax]);colorbar;
figure;imagesc(theta);
figure;imagesc(B44);
