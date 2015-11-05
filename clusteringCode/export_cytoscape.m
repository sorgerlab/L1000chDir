function export_cytoscape(t_clustered_data, LJPdir, cos_cutoff, folderlabel)
% export_cytoscape(t_nodes, LJPdir, cos_cutoff, folderlabel)

if ~exist('cos_cutoff', 'var') || isempty(cos_cutoff)
    cos_cutoff = -.02;
end
if ~exist('folderlabel', 'var') || isempty(folderlabel) 
    folderlabel = 'cytoscape/alldata';
end

Adjacency = pdist(LJPdir,'cosine');

if cos_cutoff < 0
    cos_cutoff = quantile(Adjacency, -cos_cutoff);
end

%%
AdjVals = squareform(Adjacency).*squareform(Adjacency<cos_cutoff);
[idxR,idxC, vals] = find(triu(AdjVals));

vals = 1.5 - 2*min((max(vals,.25)),.75);

%%
if ~exist(folderlabel,'dir')
    mkdir(folderlabel)
end

table2tsv(t_clustered_data(:, {'uniqueID' 'CellLine' 'DrugName' 'Conc' 'Time' ...
    'pvalue' 'pathway_role' 'cellular_function' 'DrugClass' 'Cidx' 'Cdistance' 'GRvalue'}), ...
    fullfile(folderlabel, 'nodes.tsv'))

fprintf('\nNetwork exported in %s, %i nodes, %i edges\n\n', folderlabel, ...
    sum(any(AdjVals>0)), length(idxR));

c = strcat(['uniqueID'; cellstr(t_clustered_data.uniqueID(idxR))],'\t', ...
    ['cosinedist'; num2cellstr(vals,'%.3f')],'\t', ...
    ['uniqueID1';cellstr(t_clustered_data.uniqueID(idxC))], '\n');
f = fopen(fullfile(folderlabel, 'edges.tsv'), 'w');
fprintf(f,[c{:}]);
pause(.2)
fclose(f);

f = fopen(fullfile(folderlabel, 'comments.txt'), 'w');
fprintf(f,'\nNetwork %s exported with %i nodes and %i edges\nCutoff = %.4f\n', folderlabel, ...
    sum(any(AdjVals>0)), length(idxR), cos_cutoff);
fclose(f);
