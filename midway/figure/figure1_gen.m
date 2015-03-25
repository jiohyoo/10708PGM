

acc1 = [0.4570    0.5790    0.6570    0.7300    0.7730    0.8010    0.8090    0.8300];
acc2 = [ 0.9970    0.9990    1.0000    1.0000    0.9980    1.0000    1.0000    1.0000];
ns = 250:250:2000;
plot(ns,acc1, '-ob', 'LineWidth', 1, 'MarkerSize', 7)
hold on
plot(ns, acc2, '--*r',  'LineWidth', 1, 'MarkerSize', 7)

ylim([0.4, 1.05])
	h_legend = legend('MG', 'SCGGM')

	set(h_legend,'FontSize',16);
	set(gca, 'FontSize', 13)

	xlabel('Number of training samples', 'FontSize', 16)
	ylabel('Edge recovery rate', 'FontSize', 16)


