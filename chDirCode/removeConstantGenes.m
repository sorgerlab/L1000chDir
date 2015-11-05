function [ mat,constantGenesIdx ] = removeConstantGenes( mat,thres )
%remove genes with significant low variance using the outlier function
% 
if nargin < 2
    % empirical threshold subject to adjustment for individual dataset.
    thres = 0.01;
end
logStd = log(std(mat,0,2));
totalIdx = 1:numel(logStd);

infIdx = find(isinf(logStd))';
nonInfIdx = setxor(infIdx,totalIdx);

logStd(infIdx) = [];
outlierIdx = outlier(logStd,thres);
constantGenesMetaIdx = logStd(outlierIdx) < mean(logStd);
constantGenesIdx = outlierIdx(constantGenesMetaIdx);

constantGenesIdx = [infIdx nonInfIdx(constantGenesIdx)];

mat(constantGenesIdx,:) = [];

end

