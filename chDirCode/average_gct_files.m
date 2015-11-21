function [t_samples, diffGenes, geneNames, AlldiffGenes, AllgeneNames] = average_gct_files(t_conditions, t_files, folder, VariableNames)
%
% [t_samples, diffGenes, geneNames, AlldiffGenes, AllgeneNames] = ...
%           average_gct_files(t_conditions, t_files, folder, VariableNames)
%
%	calculation of the average gene expression across conditions and plates
%
%   t_conditions:   list of unique conditions (table with set of keys)
%   t_files:        list of QNORM files (table with file names 
%                       and condition keys). Need field 'filename' and keys
%                       matching t_conditions
%   folder:         folder where the QNORM files are stored
%   VariableNames:  variables in the gct file (cdesc) to be saved and
%                       converted
%
%   t_samples:        table with conditions and statistics
%   diffGene:          characteristic directions corresponding to t_samples
%   geneNames:     genes name and order in L1000 files
%

t_samples = table;
AlldiffGenes = [];
AllgeneNames = [];



% loop throught the conditions
for i=1:height(t_conditions)
    fprintf('\n\n\t --> %i out of %i\n\n\n', i, height(t_conditions));
    
    % select the plates with at least 90 conditions measures
    t_plates = t_files( eqtable(t_conditions(i,:), t_files) & t_files.Ncond>90, :);
    
    % all replicates of the selected condition
    plates_diffGenes = [];   
    t_idx = table;
    for j=1:height(t_plates)
        gctdata = parse_gct([folder t_plates.filename{j}]);
        
        if isempty(AllgeneNames)
            AllgeneNames = gctdata.rdesc(:,7);
        else
            assert(all(strcmp(AllgeneNames, gctdata.rdesc(:,7))));
        end
        
        plates_diffGenes = [plates_diffGenes;  gctdata.mat'];
        gct_labels = cell2table(gctdata.cdesc, 'variablenames', gctdata.chd');
        t_idx = [t_idx;
            gct_labels(:,[VariableNames(:,1)' 'pert_type']) ];
    end
    
    if isvariable(t_idx, 'rna_plate')
        t_idx.rna_plate = regexprep(t_idx.rna_plate, '_X[0-9]', '');
    end
    
    t_idx.Properties.VariableNames(VariableNames(:,1)) = VariableNames(:,2);
    t_trt = group_counts(t_idx, 1:width(t_idx));
    
    %% merge the replicates of the selected conditions
    merged_diffGenes = NaN(height(t_trt),size(plates_diffGenes,2));
    for j=1:height(t_trt)
        idx = eqtable(t_trt(j,:), t_idx);
        assert(t_trt.counts(j) == sum(idx))
        
        merged_diffGenes(j,:) = mean(plates_diffGenes(idx,:),1);        
    end
    
    t_samples = [t_samples; t_trt];
    AlldiffGenes = [AlldiffGenes; merged_diffGenes];
    
end

%% merge the duplicate genes and remove -666 labels
geneNames = setdiff(unique(AllgeneNames),'-666');

diffGenes = NaN(height(t_samples), length(geneNames));
for j=1:length(geneNames)
    idx = strcmp(geneNames{j},AllgeneNames);
    diffGenes(:,j) = mean(AlldiffGenes(:,idx),2);
end

t_samples = TableToCategorical(t_samples);
