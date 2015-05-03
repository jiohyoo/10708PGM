Y = Y_te;
X = X_te;


task_num = length(X);
    mse = 0;
    
    total_sample = 0;
    for t = 1: task_num
        y_pred(:,t) = X{t} * W(:, t);
        mse = mse + sum((y_pred(:,t) - Y{t}) .^ 2);
        total_sample = total_sample + length(y_pred(:,t));
        Y_all(:,t) = Y{t};
    end
    mse = mse / total_sample;
    
    cmin = min(min(Y_all));
    cmax = max(max(Y_all));
    
    figure;imagesc(Y_all);caxis([cmin cmax]);colorbar;
     figure;imagesc(y_pred);caxis([cmin cmax]);colorbar;
    