% clear all
close all

rng('shuffle');

%%
load ResOut


%% evaluate prediction error
fprintf('---------------------------------------\n');
[prederr] = PGM_predict(opt.theta, opt.alpha1, opt.beta, opt.betad, ToyData.Y_te, ToyData.D_te);
fprintf('PGM_model (conditional) - lambda: %g\n, prediction error: %g\n', opt.lambda, prederr);

fprintf('---------------------------------------\n');
[prederr_ori] = PGM_predict(opt_ori.theta, opt_ori.alpha1, opt_ori.beta, opt_ori.betad, ToyData.Y_te, ToyData.D_te);
fprintf('PGM_model (joint) - lambda: %g\n, prediction error: %g\n', opt_ori.lambda, prederr_ori(rep));



fprintf('---------------------------------------\n');
[prederr] = PGM_predict_rnd(opt.theta, opt.alpha1, opt.beta, opt.betad, ToyData.Y_te, ToyData.D_te);
fprintf('PGM_model (conditional) - lambda: %g\n, prediction error: %g\n', opt.lambda, prederr);

fprintf('---------------------------------------\n');
[prederr_ori] = PGM_predict_rnd(opt_ori.theta, opt_ori.alpha1, opt_ori.beta, opt_ori.betad, ToyData.Y_te, ToyData.D_te);
fprintf('PGM_model (joint) - lambda: %g\n, prediction error: %g\n', opt_ori.lambda, prederr_ori(rep));

