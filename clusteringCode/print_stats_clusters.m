function print_stats_clusters(t_in, folder)


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
t_in.DrugNameTarget = categorical(strcat(cellstr(t_in.DrugName), ' (', cellstr(t_in.nominal_target), ')'));
    
for iC = 1:max(t_in.Cidx);
    
    discrete_stats_fields = {'CellLine' 'DrugNameTarget' 'Conc' 'Time' 'DrugClass' 'cellular_function'};
    cont_stats_fields = {'GRvalue' 'log10pvals'};
    add_fields = {'LJPplate' 'HMSLid'};
    
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
