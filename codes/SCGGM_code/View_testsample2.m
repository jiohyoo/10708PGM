
Theta = opt.Theta;
intercept = opt.intercept;
x_ts = X_te;
y_ts = Y_te;


K = size(Theta.yy, 2); 
N_ts = size(x_ts, 1);

Beta = scggm_indirect_SNP_overall(Theta); 

Ey_ts = x_ts * Beta  + repmat( intercept, N_ts , 1); 
if nargin > 3
	if size(y_ts, 1) ~=N_ts
		fprintf('sCGGM:error! Genotype and expression test data sample size inconsistent!\n');
		Ey_ts = []; 
		prederr = nan; 
		return; 
	end
	res  = y_ts- Ey_ts; 
	prederr = sum(sum(res.^2)) / K / N_ts;
else
	prederr = nan; 
end

cmax = max(max(Y));
cmin = min(min(Y));

figure;imagesc(y_ts);colorbar;
figure;imagesc(Ey_ts);caxis([cmin cmax]);colorbar;
figure;imagesc(Theta.yy);colorbar;
figure;imagesc(Theta.xy);colorbar;