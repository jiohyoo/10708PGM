clear all
close all

addpath(genpath('./UGM_2009')); 
addpath(genpath('./TFOCS-1.3.1'));
addpath(genpath('./pnopt-0.9-rc'));

n = 100;
p = 200;
X = randn(n,p);
y = sign( X * ( (rand(p,1) > .5) .* randn(p,1) ) + randn(n,1) );

logistic_obj = @(w) LogisticLoss(w,X,y); 
lambda = 10;
l1_pen  = prox_l1(lambda);
w0 = zeros(p,1);
pnopt_options = pnopt_optimset( 'optim', 1e-8 );

[ w, f ] = pnopt( logistic_obj, l1_pen, w0, pnopt_options );

figure
stem( w )
xlim( [ 1; p ] );