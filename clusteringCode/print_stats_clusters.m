function print_stats_clusters(t_in, folder, discrete_stats_fields, cont_stats_fields, add_fields)
% print_stats_clusters(t_in, folder, discrete_stats_fields, cont_stats_fields, add_fields))
%   one file per cluster with some statistics is saved in the folder
%   
%   t_in:       the output of clustering
%   folder:     where the reports will be saved
%   discrete_stats_fields: categorie fields to be saved (e.g. Cell
%                              line, Time, DrugClass)
%   cont_stats_fields:     continuous fields to be saved (e.g.
%                              GRvalue); quartiles are reported
%   add_fields:            additional fields reported in the list of perturbations
%



if ~exist('discrete_stats_fields', 'var') || isempty(discrete_stats_fields)
    discrete_stats_fields = {'CellLine' 'DrugNameTarget' 'Conc' 'Time' 'DrugClass' 'cellular_function'};
end
if ~exist('cont_stats_fields', 'var') || isempty(cont_stats_fields)
    cont_stats_fields = {'GRvalue' 'log10pvals'};
end
if ~exist('add_fields', 'var') || isempty(add_fields)
    add_fields = {'LJPplate' 'HMSLid'};
end

if ~exist('folder','var')
    folder = 'temp_clusters';
elseif strcmp(folder,'x')
    return
end
if ~exist(folder,'dir')
    mkdir(folder)
else
    delete([folder '/Clust_stats*txt'])
end

t_in = [t_in, table(log10(t_in.pvalue+1e-6), 'variablenames', {'log10pvals'})];
if isvariable(t_in, 'nominal_target')
    t_in.DrugNameTarget = categorical(strcat(cellstr(t_in.DrugName), ...
        ' (', cellstr(t_in.nominal_target), ')'));
else
    t_in.DrugNameTarget = t_in.DrugName;
end
    
for iC = 1:max(t_in.Cidx);
        
    file = fopen([folder '/Clust_stats_' num2str(iC) '.tsv'],'w');
    
    fprintf(file, 'Statistics for the cluster #%i\n', iC);
    fprintf(file, '---------------------------------------\n\n');
    
    fprintf(file, '#perturbations:\t%i\n\n', sum(t_in.Cidx==iC));
    
    if ~any(t_in.Cidx==iC)
        pause(.1); fclose(file);
        continue
    end
    
    %%
    for i=1:length(discrete_stats_fields)
        if iscategorical(t_in.(discrete_stats_fields{i}))
            [n,val] = hist_cat(t_in.(discrete_stats_fields{i})(t_in.Cidx==iC));
            val = cellstr(val);
        elseif isnumeric(t_in.(discrete_stats_fields{i}))
            [n,val] = hist_cat(t_in.(discrete_stats_fields{i})(t_in.Cidx==iC));
            val = cellfun2(@num2str,num2cell(val));
        end
        fprintf(file, '\n* %s distribution and count:\n', discrete_stats_fields{i});
        
        for j=1:length(val)
            fprintf(file, '\t%14s:\t%5.1f%%\t%2i\n', val{j}, 100*n(j)/sum(n), n(j));
        end
    end
    
    for i=1:length(cont_stats_fields)
        val = quantile(t_in.(cont_stats_fields{i})(t_in.Cidx==iC), [0 .25 .5 .75 1]);
        fprintf(file, '\n* %s distribution:\n', cont_stats_fields{i});
        fprintf(file, '\t[min,\t25%%\tmed\t75%%\tmax]:\n');
        fprintf(file, '\t[%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f]\n', val(:));
    end
    
    fprintf(file, '\n----------------------------------------\nall perturbations\n');
    c = table2cellstr(t_in(t_in.Cidx==iC, ...
        [discrete_stats_fields cont_stats_fields add_fields]));
    c(:,1:end-1) = strcat(c(:,1:end-1),'\t');
    for i=1:size(c,1)
        fprintf(file, [c{i,:}]);
        if i<size(c,1)
            fprintf(file, '\n');
        end
    end
    
    pause(.1);fclose(file);
    
end

fclose all;
