function [ pval] = distribution2pval(distribution,obsDist)
% calculate the pval of an experiment
% Input:
%   distribution: output of permuall.m function.
%   obsDist: the average cosine distance among chdir replicates of an 
%           experiment, calculated using: mean(pdist(expmChdirReps','cosine'))
% Output:
%   p value.
    pval = ones(numel(obsDist),1);
    for i = 1:numel(obsDist)
        pval(i) = sum(distribution<=obsDist(i))/numel(distribution);
    end

end