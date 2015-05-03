clear all
close all

load('C:\Users\HyunAh\Desktop\10708PGM\codes\HLC_data.mat');

n_sampleY = 10;
n_sampleX = 20;



Ycolsum = sum(abs(Ytrain));
[YsortV,I_Y] = sort(Ycolsum, 'descend');
Ytrain_toy = Ytrain(:,I_Y(1:n_sampleY));
figure;imagesc(Ytrain_toy);

Ytest_toy = Ytest(:,I_Y(1:n_sampleY));



Xcolsum = sum(abs(Xtrain));
[XsortV,I_X] = sort(Xcolsum, 'descend');
Xtrain_toy= Xtrain(:,I_X(1:n_sampleX));
figure;imagesc(Xtrain_toy);

Xtest_toy = Xtest(:,I_X(1:n_sampleX));


save HLC_data_toy