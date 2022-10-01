function [G,W,obj] = CDCTD(X,St,W0,StartInd,numClass,lambda)
% X is d*n centered data matrix

maxIter = 100;
numDimAfterReduction = size(W0,2);
numPts = size(X,2);
D = eye(numPts);

W = W0;
% G_init = initialClustering(X',StartInd,numClass);
G_init = full(ind2vec(StartInd',numClass))';

G = G_init;
it = 0;
obj_old = 0;
obj_new = NaN;
while ~isequal(obj_old,obj_new)&&it<maxIter
    it = it+1;
    obj_old = obj_new;    
    
    A = X'*W*W'*X;
    G_init = G;
    [G,~] = D_NCut(A,D,numClass,G_init);  
    
    Sb = X*G*(G'*G)^-1*G'*X';
    M = Sb-lambda*St;
    
    [W,~,~] = eig1(M,numDimAfterReduction);
    obj_new = trace(W'*M*W);
    obj(it) = obj_new;
end
if it == maxIter
    disp(['Warnning: the CDCTD_SP does not converge within ',num2str(maxIter),' iterations!']);
end
dbstop if error

end




