rng(0);
strings = dir('*.mat');
times = 50;

for s=1:length(strings)
    datasetName = strings(s).name;
    load(datasetName);
    numClass = length(unique(Y));
    numPts = size(X,1);
    StartInd = zeros(numPts,times);
    parfor i=1:times
        StartInd(:,i) = randsrc(numPts,1,1:numClass); 
    end

    H = eye(numPts)-ones(numPts,numPts)/numPts;
    X = X'*H;
    W0 = pca(X');
    numDimAfterReduction = min(numClass-1,size(W0,2));
    W1 = W0(:,1:numDimAfterReduction);
    
    save(datasetName,'StartInd','W1','-append');
end
    