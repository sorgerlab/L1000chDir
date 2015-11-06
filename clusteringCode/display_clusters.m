function display_clusters(t_clustered_data, Cltree, ClleafOrder, figtag, pv)
% display_clusters(t_clustered_data, Cltree, ClleafOrder, figtag, pv)
%   ## this function needs to be generalized ##
%   ## currently work for variables CellLine, Time, Conc, DrugClass,
%           GRvalue, pvalue ##
%
%   variables are the output of the clustering (function
%           L1000_chDirClustering.m)
%
%   figtag is the name of the figure
%
%   pv contains the plotting variables for the discrete variables:
%   - CellLine          (list of n names)
%   - TimePoint         (list of n names)
%   - Conc              (list of n names)
%   - DrugClass         (list of n names)
%   - CellLineColors    (nx3 color)
%   - TimePointColors   (nx3 color)
%   - ConcColors        (nx3 color)
%   - DrugClassColors   (nx3 color)
%

Cidx = t_clustered_data.Cidx;
Nclust = max(Cidx)-1;


%% distributions

DrugCellfraction = zeros(max(Cidx), length(pv.CellLine));
DrugTPfraction = zeros(max(Cidx), length(pv.TimePoint));
DrugConcfraction = zeros(max(Cidx), length(pv.Conc));
DrugCellFctfraction = zeros(max(Cidx), length(pv.DrugClass));
Clustercount = zeros(max(Cidx),1);

for i=1:max(Cidx)
    idx = t_clustered_data.Cidx==i;
    Clustercount(i) = sum(idx);
    DrugCellfraction(i,:) = hist_cat(t_clustered_data.CellLine(idx), pv.CellLine);
    DrugTPfraction(i,:) = hist(t_clustered_data.Time(idx), pv.TimePoint);
    DrugConcfraction(i,:) = hist(t_clustered_data.Conc(idx), pv.Conc);
    DrugCellFctfraction(i,:) = hist_cat(t_clustered_data.DrugClass(idx), pv.DrugClass);
end
DrugCellfraction = NormSumUnit(DrugCellfraction,2);
DrugTPfraction = NormSumUnit(DrugTPfraction,2);
DrugConcfraction = NormSumUnit(DrugConcfraction,2);
DrugCellFctfraction = NormSumUnit(DrugCellFctfraction,2);



%%


colormap(hot(12));

xpos = .09;
xwidth = .71;
xlpos = .81;
barw = .7;

ypos = [
    .95 .04  % dendo
    .77 .18  % Cell Line
    .56 .18  % drug
    .255  .1 % GRvalue
    .12  .1  % pvalue
    .465 .05 % TimePoint
    .38 .05  % Conc
    .04 .05];% cnt

get_newfigure(50, [50 50 750 950], [ figtag '.pdf'] );


%

get_newaxes([xpos ypos(1,1) xwidth ypos(1,2)])
h = dendrogram(Cltree,Nclust,'Reorder',ClleafOrder,'orientation','top');

CLpos = 1:(Nclust+1);

for i=1:length(h);set(h(i),'color','k'),end
set(gca,'xtick',[],'ytick',[],'Visible','off')
xlim([.5 Nclust+1.5])


a = get_newaxes([xpos ypos(2,1) xwidth ypos(2,2)],1);
h = bar(CLpos, DrugCellfraction, barw, 'stack');
hl = legend(h(end:-1:1), cellstr(pv.CellLine(end:-1:1)),'location','eastoutside');
posl = get(hl,'position');
set(hl,'position',[xlpos posl(2:4)])
for i=1:length(h)
    set(h(i),'facecolor',pv.CellLineColors(i,:));
end
ylim([0 1.01])


a(2) = get_newaxes([xpos ypos(3,1) xwidth ypos(3,2)],1);
h = bar(CLpos, DrugCellFctfraction, barw, 'stack');
hl = legend(h, cellstr(pv.DrugClass),'location','eastoutside','interpreter','none');
posl = get(hl,'position');
set(hl,'position',[xlpos posl(2:4)])
cnt = 0;
cols = [3 6];
ylim([0 1.01])

for i=1:length(h)
    set(h(i),'facecolor',pv.DrugClassColors(strcmp(get(h(i),'DisplayName'), pv.DrugClass),:));
end

a(3) = get_newaxes([xpos ypos(6,1) xwidth ypos(6,2)],1);
h = bar(CLpos, DrugTPfraction, barw, 'stack');
for i=1:length(h)
    set(h(i),'facecolor',(i-1)*[1 1 1]/(length(h)-1));
end
hl = legend(h, strcat(cellfun2(@num2str,num2cell(pv.TimePoint)),'h'),'location','eastoutside');
posl = get(hl,'position');
set(hl,'position',[xlpos posl(2:4)])
ylim([0 1.02])

a(4) = get_newaxes([xpos ypos(7,1) xwidth ypos(7,2)],1);
h = bar(CLpos, DrugConcfraction, barw, 'stack');
for i=1:length(h)
    set(h(i),'facecolor',pv.ConcColors(i,:));
end
hl = legend(h, strcat(cellfun2(@num2str,num2cell(pv.Conc)), '\muM'),'location','eastoutside');
posl = get(hl,'position');
set(hl,'position',[xlpos posl(2:4)])
ylim([0 1.02])

a(5) = get_newaxes([xpos ypos(4,1) xwidth ypos(4,2)],1);
plot([0 max(Cidx)+1], [1 1], '-r')
plot([0 max(Cidx)+1], [0 0], '-r')
for i=1:max(Cidx)
    plot_vbox(CLpos(i), max(-.5,t_clustered_data.GRvalue(t_clustered_data.Cidx==i)));
end
ylim([-.5 1.2])

a(6) = get_newaxes([xpos ypos(5,1) xwidth ypos(5,2)],1);
plot([0 max(Cidx)+1], -log10([.05 .05]), '-r')
plot([0 max(Cidx)+1], [3 3], '-r')
for i=1:max(Cidx)
    plot_vbox(CLpos(i), -log10(t_clustered_data.pvalue(t_clustered_data.Cidx==i)+1e-6));
end
ylim([min(-log10(t_clustered_data.pvalue))-.2 6])



a(7) = get_newaxes([xpos ypos(8,1) xwidth ypos(8,2)],1);
h = bar(CLpos, Clustercount, barw,'k');
ylim([0 max(Clustercount(1:end-1))*1.2])

for i=1:length(a)
    xlim(a(i), [.5 max(Cidx)+.5])
    set(a(i),'xtick',[])
end
set(a(7),'xtick',1:Nclust)


