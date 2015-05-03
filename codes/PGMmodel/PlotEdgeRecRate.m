clear all
close all

load SCGGM_EdgeRecRate
load EdgeRecRes Accuracy

figure;hold on;
plot([0:250:2000], [0 Accuracy'], '-rx');
plot([0:250:2000], [0 SCGGM_EdgeRecRate], '-b*');
title('Edge recovery rate');
xlabel('Number of samples');
ylabel('Edge recovery rate');
legend('MG (proposed)','SCGGM','Location', 'SouthEast');