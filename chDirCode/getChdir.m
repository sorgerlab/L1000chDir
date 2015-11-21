function [ res, expmIdx ] = getChdir( arr )
%Calculate chdir for each experiment replicate.
%   Detailed explanation goes here
ctrlIdx = cellfun(@(x)strcmp(x.pert_type,'ctl_vehicle'), arr, 'UniformOutput',false);
ctrlIdx = [ctrlIdx{:}];
expmIdx = ~ctrlIdx;

ctrlArr = arr(ctrlIdx);
ctrlMat = cellfun(@(x)x.data', ctrlArr, 'UniformOutput',false);
ctrlMat = [ctrlMat{:}];
expmArr = arr(expmIdx);

% remove outlier control replicates
ctrlOutlierIdx = pcaOutliers(ctrlMat);
ctrlMat(:,ctrlOutlierIdx) = [];
disp('Number of removed outliers in Control:');
disp(numel(ctrlOutlierIdx));

parfor i = 1:numel(expmArr)
    expmVector = expmArr{i}.data';
     %filter out constant genes.LDA would fail if such genes were keeped
    [~,constantGenesIdx] = removeConstantGenes([ctrlMat expmVector]);
    perCtrl = ctrlMat;
    perCtrl(constantGenesIdx,:) = [];
    expmVector(constantGenesIdx,:) = [];
    
    unitV = eval_chDir(perCtrl,expmVector);
    
    %insert filtered constant genes back and assign zeros to them.
    unitV = insertRemovedGenes(unitV,constantGenesIdx);
    perCtrl = insertRemovedGenes(perCtrl,constantGenesIdx);
    expmVector = insertRemovedGenes(expmVector,constantGenesIdx);
    
    expmArr{i}.chdir = unitV;
    
    %calculate projected distance
    expmArr{i}.pMetric = projectedMetric(perCtrl,expmVector,unitV);
end

res.replicateChdirs = expmArr;

    
end

