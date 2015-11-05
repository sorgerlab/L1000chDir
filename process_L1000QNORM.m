addpath('../LJP_library/chDirCode');
addpath('../LJP_library');
addpath('m:/GitHub/L1ktools/Matlab/Lib/')
addpath('d:/GitHub/L1ktools/Matlab/Lib/')

folder = './allQNORM/';
files = dir([folder 'L*gct']);
files = {files(:).name}';

%%
% get all the files and information about them
t_files = table(files, cellfun2(@(x) x{1}{1}, regexp(files, '(LJP00[5-6])_', 'tokens')), ...
    cellfun2(@(x) x{1}{1}, regexp(files, 'LJP00[5-6]_(\w+)_[0-9]*H_', 'tokens')), ...
    cell2mat(cellfun2(@(x) str2num(x{1}{1}), regexp(files, '_([0-9]*)H_', 'tokens'))), ...
    cell2mat(cellfun2(@(x) str2num(x{1}{1}), regexp(files, '_X(\d)_', 'tokens'))), ...
    cell2mat(cellfun2(@(x) str2num(x{1}{1}), regexp(files, '_n(\d*)x', 'tokens'))), ...
    'variablenames', {'filename' 'LJPplate' 'CellLine' 'Time' 'replicateCount' 'Ncond'});

VariableNames = {
    'cell_id' 'CellLine'
    'rna_plate' 'LJPplate'
    'pert_iname' 'DrugName'
    'pert_id' 'BRDid'
    'x_hmsl_id' 'HMSLid'
    'pert_dose' 'Conc'
    'pert_time' 'Time'};

% get all different conditions (cell line, LJPplate, time point)
t_conditions = unique(t_files(:,{'CellLine' 'LJPplate' 'Time'}));

t_chDir_LJP56 = table;
all_chDir_LJP56 = [];
L1000genes = [];

%%%%%%%%% calculation of the chDir across conditions and plates
%   This replace the function 'main' and explicit 'mergeChdirReplicates'


% loop throught the conditions
for i=1:height(t_conditions)
    fprintf('\n\n\t --> %i out of %i\n\n\n', i, height(t_conditions));
    
    t_plates = t_files( eqtable(t_conditions(i,:), t_files(:,{'CellLine' 'LJPplate' 'Time'})) ...
        & t_files.Ncond>100, :);
    
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
            gct_labels(expmIdx,{'cell_id' 'rna_plate' 'pert_type' ...
            'pert_iname' 'pert_id' 'x_hmsl_id' 'pert_dose' 'pert_time'}) ...
            table(cell2mat(cellfun2(@(x) x.pMetric, plateRes.replicateChdirs)), ...
            j*ones(sum(expmIdx),1),'variablenames',{'pMetric' 'plate_id'})];
    end

    t_idx.rna_plate = regexpcelltokens(t_idx.rna_plate, '^(LJP00[56])_');
    
    t_idx.Properties.VariableNames(VariableNames(:,1)) = VariableNames(:,2);
    t_trt = group_counts(t_idx, 1:(width(t_idx)-2));
    
    % merge the replicates of the selected condition
    
    % null distribution of cosine distances
    cosDistByRepCount = permuByRepCount(plates_chDir',unique(t_trt.counts));
    
    merged_chDir = NaN(height(t_trt),size(plates_chDir,2));
    pvalue = NaN(height(t_trt),1);
    ACD = pvalue;
    meanPmetric = pvalue;
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
    
    uniqueID = TableToString(t_trt(:,{'LJPplate' 'CellLine' 'Time' 'BRDid' 'Conc'}),0);
    uniqueID = table2cell(uniqueID);
    uniqueID(:,1) = strcat(uniqueID(:,1) ,'_');
    uniqueID(:,2) = strcat(uniqueID(:,2) ,'_');
    uniqueID(:,3) = strcat(uniqueID(:,3) ,'H:');
    uniqueID(:,4) = strcat(uniqueID(:,4) ,':');
    uniqueID = cellstr2str(uniqueID,'');
    
    t_chDir_LJP56 = [t_chDir_LJP56
        t_trt(:,1:(end-1)) table(pvalue, ACD, meanPmetric, uniqueID)];
    all_chDir_LJP56 = [all_chDir_LJP56; merged_chDir];
    
end
t_chDir_LJP56 = TableToCategorical(t_chDir_LJP56);

save compiled_LJP_chDir.mat L1000genes all_chDir_LJP56 t_chDir_LJP56
