function [newCidx, confMx] = align_cluster(Cidx, plotting, fixedC2, keepLastIdx)
% [newCidx, confMx] = align_cluster(Cidx, plotting, fixedC2, keepLastIdx)

if ~exist('fixedC2','var')
    fixedC2 = false;
end
if ~exist('plotting','var')
    plotting = false;
end
if ~exist('keepLastIdx','var')
    keepLastIdx = 1;
elseif ~keepLastIdx
    keepLastIdx = 0;
else
    keepLastIdx = 1;
end

newCidx = Cidx;
Nclust = max(Cidx(:))-keepLastIdx;

idx2 = cell(1,Nclust+keepLastIdx);
for j=1:Nclust+keepLastIdx
    idx2{j} = (newCidx(:,2)==j);
end
confMx = NaN(Nclust+keepLastIdx);
for i=1:Nclust+keepLastIdx
    idx1 = (newCidx(:,1)==i);
    for j=1:Nclust+keepLastIdx
        confMx(i,j) = sum( idx1 & idx2{j} )/sum( idx1 | idx2{j} );
    end
end
confMx(isnan(confMx))=1;
%% ordering the columns based on the maximum of the confusion matrix
if fixedC2
    order = 1:Nclust;
else
    [~,order] = sort(max(confMx(1:Nclust,1:Nclust)),'descend');
end
col_order = order;

%% first permurtation based on the maximum of the confusion matrix
row_order = NaN(Nclust,1);
for i=1:Nclust
    [~,row_order(i)] = max( confMx(1:Nclust, order(i)) -2*(ismember([1:Nclust]',row_order)) );
end
assert(length(unique(row_order))==Nclust);

tempconfMx = confMx(row_order, col_order);

%% try swapping pairs to imporve alignment  (independently of keepLastIdx)
col_idx = find(diag(tempconfMx)<.5);
cnt = 0;
% clf
% subplot(121)
% imagesc(tempconfMx)
SwapOccured = true;
while ~isempty(col_idx) && cnt<Nclust && SwapOccured
    i = 1; SwapOccured = false;
    while i<=length(col_idx)
        
        row_idx = find(tempconfMx(:,col_idx(i))>.3);
%         disp([col_idx(i) -1 row_idx'])
        [~,order] = sort(tempconfMx(row_idx,col_idx(i)));
        row_idx = row_idx(order);
        trConfMx = trace(tempconfMx);
        for j=1:length(row_idx)
            swap_rows = [1:Nclust]';
            swap_rows(col_idx(i)) = row_idx(j);
            swap_rows(row_idx(j)) = col_idx(i);
            if trace(tempconfMx(swap_rows,:))>(trConfMx+.1)
                temp = row_order(col_idx(i));
                row_order(col_idx(i)) = row_order(row_idx(j));
                row_order(row_idx(j)) = temp;
                tempconfMx = tempconfMx(swap_rows,:);
                SwapOccured = true;
                
%                 fprintf('Swapping %i and %i\n', col_idx(i), row_idx(j));   
%                 subplot(122)
%                 imagesc(tempconfMx)
%                 pause
%                 subplot(121)
%                 imagesc(tempconfMx)
                
                break
            end
        end
        i = i+1;
    end
    cnt = cnt+1;
end


for i=1:Nclust
    if ~fixedC2
        newCidx( Cidx(:,2)==col_order(i) ,2) = i;
    end
    newCidx( Cidx(:,1)==row_order(i) ,1) = i;
end

confMx = NaN(Nclust+keepLastIdx);
idx2 = cell(1,Nclust+keepLastIdx);
for j=1:Nclust+keepLastIdx
    idx2{j} = (newCidx(:,2)==j);
end
for i=1:Nclust+keepLastIdx
    idx1 = (newCidx(:,1)==i);
    for j=1:Nclust+keepLastIdx
        confMx(i,j) = sum( idx1 & idx2{j} )/ sum( idx1 | idx2{j} );     
    end
end
confMx(isnan(confMx))=1;

assert(all(all(confMx(1:end-keepLastIdx,1:end-keepLastIdx)==tempconfMx)))


if plotting
    cnt1 = NaN(Nclust+keepLastIdx,1);
    cnt2 = cnt1';
    for i=1:Nclust+keepLastIdx
        cnt1(i) = sum(newCidx(:,1)==i);
        cnt2(i) = sum(newCidx(:,2)==i);
    end
    
    imagesc([[cnt1/max(cnt1(1:end-keepLastIdx));0] [confMx; cnt2/max(cnt2(1:end-keepLastIdx))]],[0 1])
    title(trace(confMx)/length(confMx))
    colormap gray
end