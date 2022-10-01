function [Y,ob] = D_NCut(A,D,numClass,Yinit)
% solve Y=arg max( trace( (Y'*D*Y)^(-1)*(Y'*A*Y) ) )
% Y is n*c cluster indicator matrix

[MAinit,MDinit] = deal(zeros(numClass,1));
for i=1:numClass
    MAinit(i) = Yinit(:,i)'*A*Yinit(:,i);
    MDinit(i) = Yinit(:,i)'*D*Yinit(:,i);
end
iter = 2;
obj(1) = +inf;
obj(2) = sum(MAinit./MDinit);
Y = Yinit;MA = MAinit;MD = MDinit;
while obj(iter)~=obj(iter-1)
    iter = iter+1;
    [MA,MD,Y] = coordinateDescent(Y,A,D,MA,MD);
    obj(iter) = sum(MA./MD);
end
iter = iter-2;
obj(1) = [];
ob = obj(end);
end

function [MA,MD,Y] = coordinateDescent(Y,A,D,MA,MD)
numPts = length(Y);
for i=1:numPts
    m = find(Y(i,:)==1);
    numPts_eachClass = sum(Y);
    Y0 = Y;
    if numPts_eachClass(m)>1
        Y0(i,m) = 0;
    else
        continue;
    end
    numClass = size(Y,2);
    [MA0,MD0,MAk,MDk] = deal(zeros(numClass,1));
    for k=1:numClass
        if k==m
            MAk(k) = MA(k); MDk(k) = MD(k);
            MA0(k) = MA(k)-2*Y0(:,k)'*A(:,i)-A(i,i); MD0(k) = MD(k)-D(i,i);            
        else
            MA0(k) = MA(k); MD0(k) = MD(k);
            MAk(k) = MA(k)+2*Y(:,k)'*A(:,i)+A(i,i); MDk(k) = MD(k)+D(i,i); 
        end
    end
    delta = MAk./MDk-MA0./MD0;
    [~,p] = max(delta);
    Y = Y0;
    Y(i,p) = 1;
    if p~=m
        MA(m) = MA0(m); MD(m) = MD0(m);
        MA(p) = MAk(p); MD(p) = MDk(p);
    end
end
end

