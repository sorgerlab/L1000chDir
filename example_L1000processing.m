
addpath('./chDirCode/', './chDirCode/l1ktools_lib/')

%% get all the files and information about them
folder = './allQNORM/';
files = dir([folder 'L*gct']);
files = {files(:).name}';

% fill up the table with the information about each file (based on LJP56 file names)
t_files = table(files, ...
    cellfun2(@(x) x{1}{1}, regexp(files, '(LJP00[5-6])_', 'tokens')), ...
    cellfun2(@(x) x{1}{1}, regexp(files, 'LJP00[5-6]_(\w+)_[0-9]*H_', 'tokens')), ...
    cell2mat(cellfun2(@(x) str2num(x{1}{1}), regexp(files, '_([0-9]*)H_', 'tokens'))), ...
    cell2mat(cellfun2(@(x) str2num(x{1}{1}), regexp(files, '_X(\d)_', 'tokens'))), ...
    cell2mat(cellfun2(@(x) str2num(x{1}{1}), regexp(files, '_n(\d*)x', 'tokens'))), ...
    'variablenames', {'filename' 'LJPplate' 'CellLine' 'Time' 'replicateCount' 'Ncond'});

% get all different conditions (cell line, LJPplate, time point) based on
% the table with the file names
t_conditions = unique(t_files(:,{'CellLine' 'LJPplate' 'Time'}));

% matching the variables in the Broad L1000 .gct file to desired variable names
VariableNames = {
    'cell_id' 'CellLine'
    'rna_plate' 'LJPplate'
    'pert_iname' 'DrugName'
    'pert_id' 'BRDid'
    'x_hmsl_id' 'HMSLid'
    'pert_dose' 'Conc'
    'pert_time' 'Time'};

%% run the algorithm for the chDir
[t_chDir, chDir, L1000genes] = process_L1000QNORM(t_conditions, t_files, folder, VariableNames);

%% and save the data
save compiled_chDir.mat t_chDir chDir L1000genes
