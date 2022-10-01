close all
clear
% clc
disp('Method: CDCTD');
addpath('CDCTD');
path = 'data\';
strings = dir([path,'*.mat']);
%%%%%%%%%%%%%%%%%% result paths %%%%%%%%%%%%%%%%%%%%
d = 'results_CDCTD';
if ~exist(d,'dir')
    mkdir(d);
end

times = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:length(strings)
    datasetName = strings(s).name;
    load([path,datasetName]);
    disp(['',datasetName,'']);
    numClass = length(unique(Y));
    
    numPts = size(X,1);
    H = eye(numPts)-ones(numPts,numPts)/numPts;
    X = X'*H;
    St = X*X';
    
    lambda = 0:.1:1;
    numLambda = length(lambda);
    max_meanACC = 0;
    data = zeros(numLambda,7);
    for l=1:numLambda
        lbd = lambda(l);        
        result = zeros(times,3);
        for i=1:times
            %%%%%%%%%%%% CDCTD_SP %%%%%%%%%%%%%%%
            [G,W,~] = CDCTD(X,St,W1,StartInd(:,i),numClass,lbd);
            [Ypre,~] = vec2ind(G');
            Ypre = Ypre';
            %%%%%%%%%%%%%%%%% Measure %%%%%%%%%%%%%%%%%%%
            result(i,:) = ClusteringMeasure(Y, Ypre); % result is 1*3 vector containing Acc, NMI, purity.
        end
        mean_result = mean(result,1);
        std_result = std(result,0,1);
        data(l,:) = [lbd,mean_result,std_result];
        if data(l,2)>max_meanACC
            max_meanACC = data(l,2);
            bestLambda = lbd;
            bestResult = result;
            bestYpre = Ypre;
        end
    end
    %%%%%%%%%%%%%%%%%%% Save %%%%%%%%%%%%%%%%%%%%
    save([d,'\',datasetName],'bestLambda','bestYpre','bestResult','data');
    
end

