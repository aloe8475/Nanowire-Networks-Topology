%% Graph Theory Explore 
% This function plots graph theory parameters overlayed on graph view of
% currents for the chosen Sim at the given timestamp (IndexTime)

function [f6, f7, f8, f9, f10, f11, Explore]= graph_theory_explore(Sim,G,Adj,IndexTime,threshold_network, Explore,Graph, highlightElec, new_electrodes)
%% Participant Coefficients
if threshold_network=='t'
    return
    fprintf('Error in graph_view_threshold - you should not be seeing this');
end 
f6=figure;
currAx=gca;
p3=plot(currAx,G);
% set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p3.NodeColor='red';
p3.EdgeColor='white';
p3.NodeLabel={};

%Plot Currents
p3.MarkerSize=1.5;
p3.LineWidth=1.5;

    currs=abs(Sim.Data.Currents{IndexTime});

[j,i,~]=find(tril(Adj));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end
clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
p3.EdgeCData=cc;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);
hold on

% Plot Participant Coefficients

    p3.NodeCData=Graph.Ci;
    p_ranks=Graph.P;

edges2 = linspace(min(p_ranks),max(p_ranks),7);
bins2 = discretize(p_ranks,edges2);
p3.MarkerSize=bins2;

%Label Source
labelnode(p3,highlightElec,[new_electrodes(:).Name]); %need to make this automated.

%Need to figure out how to change colormap for Nodes seperately.


    text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P)) ' (large dot)']);

title(['Participant Coefficient Analysis Timestamp ' num2str(IndexTime)]);

%% Modular z-Score:
f7=figure;
currAx=gca;
p4=plot(currAx,G);
% set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p4.NodeColor='red';
p4.EdgeColor='white';
p4.NodeLabel={};

%Plot Currents
p4.MarkerSize=1.5;
p4.LineWidth=1.5;

    currs=abs(Sim.Data.Currents{IndexTime});

[j,i,~]=find(tril(Adj));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end
clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
p4.EdgeCData=cc;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);
hold on

    p4.NodeCData=Graph.Ci;
    mod_ranks=Graph.MZ;

edges3 = linspace(min(mod_ranks),max(mod_ranks),7);
bins3 = discretize(mod_ranks,edges3);
p4.MarkerSize=bins3;

%Label Source
labelnode(p4,highlightElec,[new_electrodes(:).Name]); %need to make this automated.


    text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ)) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);

title(['Within Module Degree z-Score Analysis Timestamp ' num2str(IndexTime)]);

%% Connectivity
f8=figure;
currAx=gca;
p5=plot(currAx,G);
% set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p5.NodeColor='red';
p5.EdgeColor='white';
p5.NodeLabel={};

%Plot Currents
p5.MarkerSize=1.5;
p5.LineWidth=1.5;

    currs=abs(Sim.Data.Currents{IndexTime});

[j,i,~]=find(tril(Adj));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end
clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
p5.EdgeCData=cc;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);
hold on


    edges = linspace(min(Graph.DEG),max(Graph.DEG),7);
    bins = discretize(Graph.DEG,edges);


p5.MarkerSize = bins;
p5.NodeColor='r';
title(['Degree Size Timestamp ' num2str(IndexTime)]);
highlight(p5,highlightElec,'NodeColor','green'); %change simulation number
labelnode(p5,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
text(-5,-6.2,'Min Degrees - 1 (small dot) | Max Degrees - 44 (large dot)');

%% Shortest Path (Distance)
d = distances(G); %calculate all the shortest path distance across all node pairs in the network.
for i = 1:length(new_electrodes)
    electrodes_cell(i)=new_electrodes(i).Name;
end

sourceIndex = find(contains(electrodes_cell,'Source')); %find index of electrodes that are source electrodes
drainIndex = find(contains(electrodes_cell,'Drain')); %find index of electrodes that are drain electrodes

sourceElec=highlightElec(sourceIndex); %change to show path from different electrodes
drainElec=highlightElec(drainIndex);
avgD=mean(d);
medianD=median(d);
stdD=std(d);

%Distribution of shortest paths:
f9=figure;
subplot(3,1,1)
h1=histogram(d);
title('Distribution of Path Distances across all Node Pairs');
xlabel('Distance');
ylabel('Frequency');
subplot(3,1,2)

h2=histogram(avgD);
title('Distribution of Mean Path Distances across all Node Pairs');
xlabel('Average Distance');
ylabel('Frequency');
subplot(3,1,3)
h3=histogram(medianD);
title('Distribution of Median Path Distances across all Node Pairs');
xlabel('Median Distance');
ylabel('Frequency');

%This needs to be fixed (automated) - what if there are more than 4 electrodes?
if length(sourceElec)>1
    source1=sourceElec(1); %choose first electrode if there are more than 1
    source2=sourceElec(2);
else
    source1=sourceElec;
end
if length(drainElec)>1
    drain1=drainElec(1);%choose first electrode if there are more than 1
    drain2=drainElec(2);
else
    drain1=drainElec;
end

%Plot all paths from current Electrode
f10=figure;
currAx=gca;
p6=plot(currAx,G);
p6.NodeLabel=d(source1,:); %label the shortest path from the source;
p6.NodeCData=d(source1,:);
highlight(p6,source1,'MarkerSize',8);
colormap(currAx,jet);%gcurrmap
colorbar(currAx);
title('Path Distances from Source Electrode');

%% Show shortest path from source to drain: %27/05/19
f11=figure;
currAx=gca;
p8=plot(currAx,G);

set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p8.NodeColor='red';
p8.EdgeColor='white';
p8.NodeLabel={};

%Plot Currents
p8.MarkerSize=1.5;
p8.LineWidth=1.5;

    currs=abs(Sim.Data.Currents{IndexTime});

[j,i,~]=find(tril(Adj));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end
clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
p8.EdgeCData=cc;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);

%plot shortest paths:
hold on
p7=plot(currAx,G);
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p7.NodeColor='green';
p7.Marker='none';
p7.EdgeColor='green';
p7.LineStyle='none';
p7.NodeLabel={};
highlight(p7,highlightElec,'NodeColor','green','Marker','o'); %change simulation number
[dist,path,pred]=graphshortestpath(Adj,source1,drain1,'Directed','false');
if length(sourceElec)>1
    [dist2,path2,pred2]=graphshortestpath(Adj,source2,drain2,'Directed','false');
    highlight(p7,path2,'EdgeColor','cyan','LineWidth',6,'LineStyle','-');
end
highlight(p7,path,'EdgeColor','cyan','LineWidth',6,'LineStyle','-');
title('Shortest Path from Sources to Drains, Overlayed on Current');


Explore.GraphView.Distances.Values=d;
Explore.GraphView.Distances.Avg=avgD;
Explore.GraphView.Distances.Std=stdD;
Explore.GraphView.Distances.Median=medianD;
Explore.GraphView.Distances.DistancesFromSource=d(sourceElec,:);
Explore.GraphView.Distances.ShortestPath=path;
Explore.Electrodes.Source=sourceElec;
Explore.Electrodes.Drain=drainElec;

end 