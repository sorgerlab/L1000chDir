function [newCidx, Ccenters, Cdistance, ClDist, Cltree, ClleafOrder] = L1000_chDirClustering(chDir, clusteringPara)
% [newCidx, Ccenters, Cdistance, ClDist, Cltree, ClleafOrder] = L1000_chDirClustering(chDir, clusteringPara)
% 
% Nrepeat = 101; 
% Cutoff = 0.55;
% pval = .05;
% opt = [1.22, 150, 2e-6, 0];
% 
% clusteringPara = struct('Nrepeat', Nrepeat, 'Cutoff', Cutoff, ...
%     'maxpval', pval, 'FuzzyOpt', opt);
%

Nclust = clusteringPara.Nclust;
Nrepeat = clusteringPara.Nrepeat;
Cutoff = clusteringPara.Cutoff;

wCidx = NaN(Nclust, size(chDir,1), Nrepeat);
oldCidx = NaN(size(chDir,1), Nrepeat);
maxCidx = oldCidx;

% iterative process
parfor i=1:Nrepeat
    s = RandStream('mt19937ar','Seed',i);
    RandStream.setGlobalStream(s);
    [~, temp, obj_fct] = fcm_cosineDist(chDir, Nclust, clusteringPara.FuzzyOpt);
    [maxCidx(:,i), oldCidx(:,i)] = max(temp);
    oldCidx(:,i) = oldCidx(:,i).*(max(temp)>=Cutoff)' ...
        + (Nclust+1)*(max(temp)<Cutoff)';
    fprintf('%.1f ', obj_fct(end));
    if mod(i,20)==0, fprintf('\n'); end
    wCidx(:,:,i) = temp;
end
fprintf('\n');

% disply the clustering
figure(1001)
clf
subplot(121)
imagesc(wCidx(:,:,1),[0 1])
colormap gray

subplot(122)
hist(max(wCidx(:,:,1)),0:.05:1)
mean(mean(max(wCidx)>=Cutoff))

%% alignment of the clusters
distcluster = NaN(Nrepeat);
for i=1:Nrepeat
    parfor j=(i+1):Nrepeat
        [~, confMx] = align_cluster(oldCidx(:,[i j]));
        distcluster(i,j) = trace(confMx);
    end
end

for i=1:Nrepeat
    for j=(i+1):Nrepeat
        distcluster(j,i) = distcluster(i,j);
    end
end
[~,Cref] = max(nansum(distcluster));
[~,Cref2] = nanmax(distcluster(:,Cref));


% display the results of the alignment
figure(1002)
clf
subplot(10,10,1)
newCidx = NaN*oldCidx;
allTraces = NaN(1,Nrepeat-1);
[newCidx(:,[1 2]), confMx] = align_cluster(oldCidx(:,[Cref Cref2]),1);
allTraces(1) = trace(confMx);
%
[~,order] = sort(distcluster(:,Cref),'descend');
for i=3:Nrepeat
    subplot(10,10,i-1)
    [temp, confMx] = align_cluster([oldCidx(:,order(i)) newCidx(:,1)],1,1);
    assert(all(temp(:,2)==newCidx(:,1)))
    newCidx(:,i) = temp(:,1);
    allTraces(i-1) = trace(confMx);
end

% find the cluster for each perturbation
Cidx = (Nclust+1)*ones(size(newCidx,1),1);
for i=1:size(newCidx,1)
    n = hist(newCidx(i,:),1:(Nclust+1));
    if max(n)>(Nrepeat/2);
        [~,Cidx(i)] = max(n);
    end
end
Cidx = max_cntvalue(newCidx');
allCidx = newCidx;


%% Evaluate the center of each cluster
Ccenters = NaN(Nclust, size(chDir,2));
Cdistance = NaN(size(chDir,1),1);
for i=1:Nclust
    Ccenters(i,:) = mean(chDir(Cidx==i,:),1);
    Ccenters(i,:) = Ccenters(i,:)/norm(Ccenters(i,:));
    
    Cdistance(Cidx==i) = 1-(chDir(Cidx==i,:)*Ccenters(i,:)');
end



%% annotation of the perturbations

Cidx = ToColumn(Cidx);


ClDist = 1-Ccenters(1:Nclust,:)*Ccenters(1:Nclust,:)';
ClDist = ClDist.*(1-eye(size(ClDist)));
Cltree = linkage(ClDist,'average');
ClDist(eye(size(ClDist))==1) = 0; 
ClDist(isnan(ClDist)) = 2;
CLsqrDist = squareform(ClDist);
ClleafOrder = optimalleaforder(Cltree,CLsqrDist);

newCidx = max(Cidx)*ones(size(Cidx));
for i=1:length(ClleafOrder)
    newCidx(ClleafOrder(i)==Cidx) = i;
end

Ccenters = Ccenters(ClleafOrder,:);
ClDist = ClDist(ClleafOrder,ClleafOrder);
