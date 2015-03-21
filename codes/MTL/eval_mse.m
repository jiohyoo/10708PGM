function mse = eval_mse (Y, X, W)
% Evaluate Mean squared error
%
    task_num = length(X);
    mse = 0;
    
    total_sample = 0;
    for t = 1: task_num
        y_pred = X{t} * W(:, t);
        mse = mse + sum((y_pred - Y{t}) .^ 2);
        total_sample = total_sample + length(y_pred);
    end
    mse = mse / total_sample;
end
