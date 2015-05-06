clear all
close all

load EdgeRecoveryRate_3state_n2000_rep1

% thre = 2% 0.85;

figure;imagesc(ToyData.maskDisCts);colorbar;
rep = 1;
for ii = 1 : size(opt,1)
    
    thetaT = opt{ii, rep}.theta';
    figure;imagesc(thetaT);
    IDX = kmeans(thetaT(:), 2);
    theta_edges = zeros(size(thetaT,1), size(thetaT,2));
%     ooo = opt{ii,rep}.theta(IDX);
%     oooo = reshape(ooo, size(opt{ii, rep}.theta,1),size(opt{ii, rep}.theta,2));
%     theta_edges = ooofo';
    IDX1 = find(IDX==1);
    IDX2 = find(IDX==2);
    theta_edges(IDX1) = 1;
    theta_edges(IDX2) = 0;
    if thetaT(IDX1(1)) > thetaT(IDX2(1))
        theta_edges(IDX1) = 1;
        theta_edges(IDX2) = 0;
    else
        theta_edges(IDX1) = 0;
        theta_edges(IDX2) =1;
    end
    figure;imagesc(theta_edges);
    
%     theta_edges = double((opt{ii, rep}.theta' > thre )); % binary
    theta_edges2 = zeros(size(theta_edges,1), size(theta_edges,2)/3);
    for jj = 1: size(theta_edges,2)/3
        theta_edges2(:,jj) = sum(theta_edges(:,(jj-1)*3+1:(jj-1)*3+3),2); % 2 states - OR op
    end
    theta_edges2 = theta_edges2 > 0;
    
    TP(ii,rep) = sum(sum(theta_edges2 .* ToyData.maskDisCts));
    FP(ii,rep) = sum(sum(theta_edges2 .* (ToyData.maskDisCts * (-1) + 1)));
    TN(ii,rep) = sum(sum((theta_edges2 * (-1) + 1) .* (ToyData.maskDisCts * (-1) + 1)));
    FN(ii,rep) = sum(sum((theta_edges2 * (-1) + 1) .* ToyData.maskDisCts));
    
    
    figure;imagesc(opt{ii, rep}.theta');colorbar;
    figure;imagesc(theta_edges);colorbar;
    figure;imagesc(theta_edges2);colorbar;
    figure;hist(opt{ii, rep}.theta(:));
%     mean(opt{ii, rep}.theta(:))
%     std(opt{ii, rep}.theta(:))
%     median(opt{ii, rep}.theta(:))

    
end



%% plot!
TP_avg = mean(TP,2);
FP_avg = mean(FP,2);
TN_avg = mean(TN,2);
FN_avg = mean(FN,2);
Accuracy = (TP_avg + TN_avg) ./ (TP_avg + FP_avg + TN_avg + FN_avg);
figure;
plot([0:250:2000], [0; Accuracy]);

