function arr = convert_gct_arr(gctdata)

%%
arr = cell(length(gctdata.cid),1);
for i=1:length(gctdata.cid)
    arr{i}.data = gctdata.mat(:,i)';
    arr{i}.pert_type = gctdata.cdesc(i,27);
end
