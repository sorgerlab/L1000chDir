function [ chdirs ] = mergeChdirReplicates( platesRes,sigIdStructs )
%Merge chdir replicates and computes the significance.
%   Detailed explanation goes here

chdirStructsAllPlates = cellfun(@(x)x.replicateChdirs',platesRes,'UniformOutput',false);
chdirStructsAllPlates = [chdirStructsAllPlates{:}];
chdirReplicatesTotal = cellfun(@(x)x.chdir, chdirStructsAllPlates,'UniformOutput',false);
chdirReplicatesTotal = [chdirReplicatesTotal{:}];

repCount = cellfun(@(x)x.replicateCount,sigIdStructs,'UniformOutput',false);
repCount = [repCount{:}];

% caculate cosine distance null distribution.
% cosDist = permuall(chdirReplicatesTotal,repCount);
cosDist = permuByRepCount(chdirReplicatesTotal,repCount);



chdirs = cell(numel(sigIdStructs),1);
for i = 1:numel(sigIdStructs)
    sigIdStruct = sigIdStructs{i};
    chdir = sigIdStruct.x0x5F_id;
    chdir.replicateCount = sigIdStruct.replicateCount;
    if sigIdStruct.replicateCount==1
        chdirReplicate = dict(chdirStructsAllPlates,@(x)strcmp(x.pert_id,sigIdStruct.x0x5F_id.pert_id)&&strcmp(x.pert_dose,sigIdStruct.x0x5F_id.pert_dose));
        plateRes = dict(platesRes,@(x)strcmp(chdirReplicate.det_plate,x.replicateChdirs{1}.det_plate));
        chdir.chdir = chdirReplicate.chdir';
        chdir.pvalue =  distribution2pvalRight(plateRes.ctrlPMetrics, chdirReplicate.pMetric);
%         chdir.zscore = (chdirReplicate.pMetric-mean(plateRes.ctrlPMetrics))/std(plateRes.ctrlPMetrics);
        chdir.ACD = chdirMeanDist;
        chdirs{i} = addFields(chdir,chdirReplicate);
    else
        chdirReplicates = dict(chdirStructsAllPlates,@(x)strcmp(x.pert_id,sigIdStruct.x0x5F_id.pert_id)&&strcmp(x.pert_dose,sigIdStruct.x0x5F_id.pert_dose));
        chdirVectors = cellfun(@(x)x.chdir,chdirReplicates,'UniformOutput',false);
        chdirVectors = [chdirVectors{:}];
        meanChdirVector = mean(chdirVectors,2);
        chdir.chdir = (meanChdirVector/norm(meanChdirVector))';
        chdirMeanDist = mean(pdist(chdirVectors','cosine'));
        eval(sprintf('cosDistByRepCount = cosDist.repCount%d;',chdir.replicateCount));
        chdir.pvalue = distribution2pval(cosDistByRepCount,chdirMeanDist);
%         chdir.pvalue = distribution2pval(cosDist,chdirMeanDist);
%         chdir.zscore = (mean(cosDist)-chdirMeanDist)/std(cosDist);
        disp(i);
        chdir.ACD = chdirMeanDist;
        chdirs{i} = addFields(chdir,chdirReplicates{1});
    end
end

end

function [chdir] = addFields(chdir,proto)
    chdir.batch = proto.batch;
    chdir.pert_id = proto.pert_id;
    chdir.pert_dose = proto.pert_dose;
    chdir.sig_id = strjoin({proto.batch,proto.pert_id,proto.pert_dose},':');
    
    chdir.pert_iname = proto.pert_iname;
    chdir.pert_type = proto.pert_type;
    chdir.pert_time = proto.pert_time;
    chdir.pert_time_unit = proto.pert_time_unit;
    chdir.pert_dose_unit = proto.pert_dose_unit;
    chdir.cell_id = proto.cell_id;
    
end

function [elems] = dict(arr,matchFunc)

   matchIdx = cellfun(matchFunc,arr,'UniformOutput',false);
   matchIdx = [matchIdx{:}];
   elems = arr(matchIdx);
   if numel(elems)==1
       elems = elems{1};
   end

end