function [ theOutliers ] = pcaOutliers( X,alpha )
%identify outliers in gene expression matrix. It is specifically designed 
%for outliers that are observable in PCA plotting.
%   In X, rows are genes and columns are samples.
if nargin < 2
    alpha = 0.01;
end
dist = pdist(X','euclidean');
distMat = squareform(dist);
medians = median(distMat);
theOutliers = outlier(medians,alpha);



end

