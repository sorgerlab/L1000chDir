function [ mat ] = insertRemovedGenes( mat,rmIdx )
%nsert a row of zeros for removed genes
%   Detailed explanation goes here
rmIdx = sort(rmIdx);
for i = 1:numel(rmIdx)
    mat = insertRow(mat,rmIdx(i)); 
end

end

function [matInserted] = insertRow(mat,idx)
    row = zeros(1,size(mat,2));
    matInserted = [mat(1:idx-1,:);row;mat(idx:end,:)];

end