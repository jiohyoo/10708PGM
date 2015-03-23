function [X_te, Y_te, D_te] = GenerateTestSample(ToyData)

n = ToyData.n_te;
p = ToyData.p;
q = ToyData.q;
thcts = ToyData.thcts;
maskDisCts = ToyData.maskDisCts;
thdiscts= ToyData.thdiscts;
L = ToyData.L;


%% generate dis var
nodePot = ToyData.nodePot;
edgePot = ToyData.edgePot;
edgeStruct = ToyData.edgeStruct;
edgeStruct.maxIter = n;
burnIn = ToyData.burnIn;
Y = UGM_Sample_Gibbs(nodePot,edgePot,edgeStruct,burnIn);
Y =Y';

Y_te = Y;



%% generate conti var

Sigma=inv(thcts);
X= zeros(n,p); 
R=chol(Sigma,'lower');
for i=1:n
    idx=Y(i,:);
    phi=zeros(p,1); % this is a funtion of idx i.e. discrete labels
    ind=find(maskDisCts);
    for ii=1:length(ind)
        [s j]=ind2sub(size(maskDisCts),ind(ii)); % location of ind (1 in maskDisCts) in coordinate
        thsj=thdiscts{s,j};
        phi(s)=phi(s)+thsj(idx(j));
    end
    mu=Sigma*phi;
    X(i,:)=(mu+R*randn(p,1))';
end

X_te = X;






D=[];
for j=1:q
    Dj=zeros(n,L(j));
    for i=1:n
        Dj(i,Y(i,j))=1;
    end
    D=[D Dj];
end

D_te = D;