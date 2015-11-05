function [ pMetrics,unitVTotal ] = ctrlPermuPmetric( ctrl,k )
%generate null distribution of projected distances by only sampling control
%replicates.
%
% Input:
%   ctrl is a matrix(978*n) of normlized 978-gene control replicates.
%   k is the prevalent experiment replicate count (if most experiments have 4 replicates, k will be 4).
% Output:
%  pMetrics is a vector of projected distances that represent the null distribution.
%  unitVTotal is an optional output that store all the computed null chdirs.
%
% p-value of an experiment could be easily computed by compare its projected 
% distance to this pMetrics, the null distriubtion of projected distances.
%
disp('Begin pMetric control permutation...')
ctrl = removeConstantGenes(ctrl);
sampleCount = size(ctrl,2);
sampleIdx = 1:sampleCount;
combinations = nchoosek(sampleIdx,k);
combinationsCount = size(combinations,1);

% Perform 10000 permutations at most.
maxi = 10000;
if combinationsCount > maxi
    randChosenIdx = randi(combinationsCount,maxi,1);
    combinations = combinations(randChosenIdx,:);
    combinationsCount = size(combinations,1);
end

pMetrics = zeros(combinationsCount,1);
unitVTotal = zeros(size(ctrl,1),combinationsCount);
parfor i = 1:combinationsCount
    combination = combinations(i,:);
    pseudoExpIdx = setxor(sampleIdx,combination);
    pseudoCtrl = ctrl(:,combination);
    pseudoExp = ctrl(:,pseudoExpIdx);
    unitV = chDir(pseudoCtrl,pseudoExp);
    pMetrics(i) =  projectedMetric(pseudoCtrl,pseudoExp,unitV);
    unitVTotal(:,i) = unitV;
end
disp('Permuation ends')
end

