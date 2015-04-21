clear all
close all

%% train 
Dir = 'C:\Users\HyunAh\Desktop\data_hlc\hlc\';
fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-xtrain.txt']);
[A,COUNT] = fscanf(fid,'%f');
fclose(fid);

fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-xtrain.txt']);
xx = fgetl(fid);
C = strsplit(xx,' ');
fclose(fid);

Ncol = size(C,2) - 1;
Nrow = size(A,1) / Ncol;

for i = 1: Nrow
    Xtrain(i,:) = A((i-1)*Ncol+1: (i-1)*Ncol+Ncol);
end

fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-ytrain.txt']);
[A4,COUNT] = fscanf(fid,'%f');
fclose(fid);
fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-ytrain.txt']);
xx4 = fgetl(fid);
C4 = strsplit(xx4,' ');
fclose(fid);
Ncol4 = size(C4,2) - 1;
Nrow4 = size(A4,1) / Ncol4;
for i = 1: Nrow4
    Ytrain(i,:) = A4((i-1)*Ncol4+1: (i-1)*Ncol4+Ncol4);
end
%% test
fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-xtest.txt']);
[A2,COUNT2] = fscanf(fid,'%f');
fclose(fid);

fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-xtest.txt']);
xx2 = fgetl(fid);
C2 = strsplit(xx2,' ');
fclose(fid);

Ncol2 = size(C2,2) - 1;
Nrow2 = size(A2,1) / Ncol2;

for i = 1: Nrow2
    Xtest(i,:) = A2((i-1)*Ncol2+1: (i-1)*Ncol2+Ncol2);
end

fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-ytest.txt']);
[A3,COUNT3] = fscanf(fid,'%f');
fclose(fid);
fid = fopen([Dir 'geno-chr6-thr0.05-K100-std-ytest.txt']);
xx3 = fgetl(fid);
C3 = strsplit(xx3,' ');
fclose(fid);
Ncol3 = size(C3,2) - 1;
Nrow3 = size(A3,1) / Ncol3;
for i = 1: Nrow3
    Ytest(i,:) = A3((i-1)*Ncol3+1: (i-1)*Ncol3+Ncol3);
end

%% processing
Xtrain_raw = Xtrain;
Xtest_raw = Xtest;

[Xtrain, Xtrain_uni, NXtrain_uni] = cvt2States(Xtrain);
[Xtest, Xtest_uni, NXtest_uni] = cvt2States(Xtest);

figure;
hist(NXtrain_uni);
figure;
hist(NXtest_uni);

idx = find(NXtest_uni == 1);
for i = 1: length(idx)
    disp(Xtest_uni{idx(i)});
end

figure;
hist(Xtrain_raw(:,847)); % 847th gene has values {-1, 0, 1}
figure;
hist(Xtest_raw(:,487)); % 487th gene has 2 states {0, 2} no {1}.
% all the others that have 2 states are comprised of two states of {0,1}
% excluding {2}.
figure;
hist(Xtest(:,487));
save HLC_data Xtrain_raw Xtest_raw Xtrain Xtest Ytest Ytrain




%% additional checks!
for i = 1: size(Xtrain_raw,2)
   nuni(i) = length(unique(Xtrain_raw(:,i)));
end

idx = find(nuni~=3);
nuni_non3 = nuni(idx);

for i = 1: length(nuni_non3)
    uni3{i} = unique(Xtrain_raw(:,idx(i)));
    if uni3{i}(2) > 1
        disp('2nd val > 1');
    else
%         disp('2nd val < 1');
    end
end


figure;
hist(Xtrain_raw(:,1));

figure;
hist(Xtrain_raw(:,900));