function [t_chDir, chDir, L1000genes] = process_L1000QNORM(t_conditions, t_files, folder, VariableNames)
%
% [t_chDir, chDir, L1000genes] = process_L1000QNORM(t_conditions, t_files, folder, VariableNames)
%
%	calculation of the chDir across conditions and plates (see example_L1000processing.m)
%
%   t_conditions:   list of unique conditions (table with set of keys)
%   t_files:        list of QNORM files (table with file names 
%                       and condition keys). Need field 'filename' and keys
%                       matching t_conditions
%   folder:         folder where the QNORM files are stored
%   VariableNames:  variables in the gct file (cdesc) to be saved and
%                       converted
%
%   t_chDir:        table with conditions and statistics
%   chDir:          characteristic directions corresponding to t_chDir
%   L1000genes:     genes name and order in L1000 files
%

t_chDir = table;
chDir = [];
L1000genes = [];



% loop throught the conditions
for i=1:height(t_conditions)
    fprintf('\n\n\t --> %i out of %i\n\n\n', i, height(t_conditions));
    
    % select the plates with at least 90 conditions measures
    t_plates = t_files( eqtable(t_conditions(i,:), t_files) & t_files.Ncond>90, :);
    
    % all replicates of the selected condition
    plates_chDir = [];   
    t_idx = table;
    for j=1:height(t_plates)
        gctdata = parse_gct([folder t_plates.filename{j}]);
        
        if isempty(L1000genes)
            L1000genes = gctdata.rdesc(:,7);
        else
            assert(all(strcmp(L1000genes, gctdata.rdesc(:,7))));
        end
        
        [plateRes, expmIdx] = getChdir(convert_gct_arr(gctdata));
        plates_chDir = [plates_chDir; 
            cell2mat(cellfun2(@(x) x.chdir', plateRes.replicateChdirs))];
        gct_labels = cell2table(gctdata.cdesc, 'variablenames', gctdata.chd');
        t_idx = [t_idx;
            gct_labels(expmIdx,[VariableNames(:,1)' 'pert_type']) ...
            table(cell2mat(cellfun2(@(x) x.pMetric, plateRes.replicateChdirs)), ...
            j*ones(sum(expmIdx),1),'variablenames',{'pMetric' 'plate_id'})];
    end
    
    if isvariable(t_idx, 'rna_plate')
        t_idx.rna_plate = regexprep(t_idx.rna_plate, '_X[0-9]', '');
    end
    
    t_idx.Properties.VariableNames(VariableNames(:,1)) = VariableNames(:,2);
    t_trt = group_counts(t_idx, 1:(width(t_idx)-2));
    
    %% merge the replicates of the selected condition
    
    % null distribution of cosine distances
    cosDistByRepCount = permuByRepCount(plates_chDir',unique(t_trt.counts));
    
    merged_chDir = NaN(height(t_trt),size(plates_chDir,2));
    pvalue = NaN(height(t_trt),1);
    ACD = pvalue;
    meanPmetric = pvalue;
    Repcout = pvalue;
    parfor j=1:height(t_trt)
        idx = eqtable(t_trt(j,:), t_idx);
        assert(t_trt.counts(j) == sum(idx))
        
        if t_trt.counts(j)==1 % only one replicate
            merged_chDir(j,:) = plates_chDir(idx,:);
            pvalue(j) = NaN;
            ACD(j) = NaN;
            meanPmetric(j) = t_idx.pMetric(idx);
        else
            meanChdirVector = mean(plates_chDir(idx,:),1);
            merged_chDir(j,:) = meanChdirVector/norm(meanChdirVector);
            chdirMeanDist = mean(pdist(plates_chDir(idx,:),'cosine'));
            pvalue(j) = distribution2pval(...
                cosDistByRepCount.(['repCount' num2str(t_trt.counts(j))]),...
                chdirMeanDist);
            ACD(j) = chdirMeanDist;
            meanPmetric(j) = mean(t_idx.pMetric(idx));
        end
        
    end
    
    t_chDir = [t_chDir
        t_trt(:,1:(end-1)) table(pvalue, ACD, meanPmetric) t_trt(:,end)];
    chDir = [chDir; merged_chDir];
    
end
t_chDir = TableToCategorical(t_chDir);