
addpath('./clusteringCode/')

%% clustering all perturbations

% parameters for clustering
clusteringPara = struct(...
    'Nrepeat', 50, ...
    'Nclust', 5, ...
    'Cutoff', .55, ...
    'maxpval', .05, ...
    'FuzzyOpt', [1.22, 150, 2e-6, 0]);

% selecting significant perturbations
Selected_chDir = chDir(t_chDir.pvalue<clusteringPara.maxpval,:);
t_selected = t_chDir(t_chDir.pvalue<clusteringPara.maxpval,:);

% clustering algorithm
[Cidx, Ccenters, Cdistance, ClDist, Cltree, ClleafOrder] = ...
    L1000_chDirClustering(Selected_chDir, clusteringPara);

% annotation of the perturbation
t_clustered_data = [t_selected  table(Cidx,Cdistance)];

%% saving the result
save(['clustered_L1000_p' num2str(clusteringPara.maxpval) '_Ncl' num2str(clusteringPara.Nclust) '.mat'], ...
    't_clustered_data', 'Selected_chDir', 'Ccenters', 'clusteringPara',  ...
     'ClDist', 'Cltree', 'ClleafOrder');
 
%% print the results for each cluster

% this should be fill up with relevant information
%   these are default output; not everything is needed (see help for
%   print_stats_clusters)

t_clustered_data.DrugClass = ...
    categorical(repmat({'-'}, height(t_clustered_data),1));
t_clustered_data.cellular_function = ...
    categorical(repmat({'-'}, height(t_clustered_data),1));
t_clustered_data.GRvalue = rand(height(t_clustered_data),1);

savefolder = ['clusters_stats_Ncl' num2str(clusteringPara.Nclust)];

print_stats_clusters(t_clustered_data, savefolder)
table2tsv(t_clustered_data, [savefolder '.tsv'])

%% displaying the results with bars

% need a structure with parameters for plotting 
%     ## this function needs to be generalized ##
%     ## currently working for the variables:
%               CellLine, TimePoint, Conc, DrugClass, GRvalue
%
pv = struct( ...
    'CellLine' , unique(t_clustered_data.CellLine), ...
    'TimePoint', unique(t_clustered_data.Time), ...
    'Conc', unique(t_clustered_data.Conc), ...
    'DrugClass', unique(t_clustered_data.DrugClass), ...
    'CellLineColors' , jet(length(unique(t_clustered_data.CellLine))), ...
    'TimePointColors', gray(length(unique(t_clustered_data.Time))), ...
    'ConcColors', hot(length(unique(t_clustered_data.Conc))), ...
    'DrugClassColors', prism(length(unique(t_clustered_data.DrugClass))));
    
display_clusters(t_clustered_data, Cltree, ClleafOrder, savefolder, pv)


