

n_iter = 10
num_samples=250:250:2000;
edge_recover_acc = zeros(size(num_samples))

%% sCGGM demo with cross-validation 
% specify the search grid of regularization parameters
lambda1_seq = [1 0.64 0.32 0.16 0.08 0.04 0.02]; 
lambda2_seq = [1 0.64 0.32 0.16 0.08 0.04 0.02]; 

% performs kcv-fold cross validation, kcv must be >= 3
kcv = 5;
% loading data
data = load('./../ToyData.mat');
X = data.Y; % SNPs
Y = data.X; % expression rate
n = data.n; % number of samples


for num_samples_idx=1:length(num_samples)
		num_sample = num_samples(num_samples_idx);
		
		% split training / test set 
		training_ratio = 0.8;
	
		num_acc = 0;
		for iter=1:n_iter
			rand_perm = randperm(n);
	
			X_tr = X(rand_perm(1:num_sample), :);
			Y_tr = Y(rand_perm(1:num_sample), :);
			
			
			%fprintf('sCGGM running ToyData with %d-fold cross-validation...\n', kcv);
			
			option.verbose = true; 
			option.maxiter = 500; 
		
			opt = scggm_cv(X_tr, Y_tr, kcv, lambda1_seq, lambda2_seq, option);
			
			num_acc = num_acc + sum(sum((opt.Theta.xy == 0) == (data.maskDisCts == 0)));
		end
		edge_recover_acc(num_samples_idx) = num_acc / (n_iter * prod(size(opt.Theta.xy)))
end
