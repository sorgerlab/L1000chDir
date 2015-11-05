function [ res ] = permuByRepCount( unitV,repCount )
% calculate empirical null distribution of chdir repliacates by 100000 permutation 
% Input:
%   unitV: (978*n) matrix, all chdir replicates of all experiments in a batch,
%                    n is the number of chdir reps of all experiments in a batch.
%   repCount: vector of length m, store replicate count of each experiment
%              in the batch. m equals to the number of experiments in the
%              batch
%
% generate a distributioni for each count >=2

permuCount = 50000;
[repCountUnique,~] = count_unique(repCount);
noDistIdx = repCountUnique == 1; % no between-rep distance if there is only one rep.
repCountUnique(noDistIdx) = [];
expmCount = size(unitV,2);

for i = 1:numel(repCountUnique)
    currentRepCount = repCountUnique(i);
    eval(sprintf('res.repCount%d = zeros(permuCount,1);',currentRepCount));
    for j = 1:permuCount
        permu = randperm(expmCount,currentRepCount);
        sample = unitV(:,permu);
        distVal = mean(pdist(sample','cosine'));
        eval(sprintf('res.repCount%d(j) = distVal;',currentRepCount));
    end
end




