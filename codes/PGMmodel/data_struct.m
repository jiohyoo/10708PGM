function [Data] = data_struct(Xtrain, Ytrain, Xtest, Ytest)

Data.p = size(Ytrain,2);
Data.q = size(Xtrain,2);
Data.n = size(Xtrain,1);
Data.L = 3*ones(size(Xtrain,2),1);

Dtr = zeros(Data.n, Data.L(1) * Data.p);
for i = 1: Data.n
    for j = 1: Data.q
        Dtr(i, (j-1)*Data.L(1)+Xtrain(i,j)+1) = 1;
    end
end

Dte = zeros(size(Xtest,1), Data.L(1) * size(Xtest,2));
for i = 1: size(Xtest,1)
    for j = 1: size(Xtest,2)
        Dte(i, (j-1)*Data.L(1)+Xtest(i,j)+1) = 1;
    end
end

Data.Y_tr = Xtrain;
Data.X_tr = Ytrain;
Data.D_tr = Dtr;
Data.Y_te = Xtest;
Data.X_te = Ytest;
Data.D_te = Dte;

Data.n_tr = size(Ytrain,1);