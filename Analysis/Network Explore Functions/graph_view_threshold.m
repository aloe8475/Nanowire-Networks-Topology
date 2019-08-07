%% Graph View
% This function plots graph parameters such as current, voltage and
% resistance for the chosen Sim at the given timestamp (IndexTime)
function [f3, f4, f5, G, Adj, Adj2, Explore,  highlightElec, new_electrodes]= graph_view_threshold(Sim,Graph,IndexTime,Explore,G, threshold_network, threshold, drain_exist,node_indices)
if threshold_network~='t'
    return
    fprintf('Error in graph_view_threshold - you should not be seeing this');
end
%Plot Graph
f3=figure;
currAx=gca;
p=plot(currAx,G); %G = binarised graph
% set(currAx,'Color',[0.35 0.35 0.35]); %change background
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p.NodeColor='red';
p.EdgeColor='white';
p.NodeLabel={};

%% Plot Currents
p.MarkerSize=1.5;
p.LineWidth=1.5;
Adj=(Sim.Data.AdjMat{IndexTime});%we need to keep a copy of the original Adj matrix (unthresholded) to find all the currents - this is the NON BINARISED Adj matrix
Adj2=Adj(threshold,threshold); %this is a different adj mat (unbinarised) compared to Graph.AdjMat (binarised by resistance)
if ~drain_exist %if no drains
    currs=log10(abs(Sim.Data.Currents{IndexTime}));
else
    currs=abs(Sim.Data.Currents{IndexTime}); %Graph.BinarisedCurrents);

end

%% What we are doing here is finding the adj matrix, and finding the edges that have current flowing through them.
%Find currents
currs=currs(threshold,threshold); %currs(node_indices);%
[j,i,~]=find(tril(Adj2));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
[j,i,~]=find(tril(Adj2));
cc2=zeros(1,length(j));

%Find edges in Adj matrix that have current in them
for k=1:length(j)
    cc2(k)=Graph.networkThreshold(i(k),j(k));
end
% remove edges in adj matrix that don't have current
cc3=cc(logical(cc2));


clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
p.EdgeCData=cc3;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);

%Highlight Electrodes:
node_indices=find(threshold==1); %find nodes with threshold == 1
for i=1:size(Sim.Electrodes.PosIndex,1)
    if ~isempty(find(node_indices==Sim.Electrodes.PosIndex(i)))
        new_electrodes(i).PosIndex=find(node_indices==Sim.Electrodes.PosIndex(i));
        new_electrodes(i).Name=Sim.Electrodes.Name(i);
    end
end

highlightElec={new_electrodes.PosIndex};
highlightElec=cell2num(highlightElec);
highlight(p,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
labelnode(p,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
title(['Current Graph View | T= ' num2str(IndexTime)]);


%Plot Graph
f4=figure;
currAx=gca;
p1=plot(currAx,G);
% set(currAx,'Color',[0.35 0.35 0.35]); %change background
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p1.NodeColor='red';
p1.EdgeColor='white';
p1.NodeLabel={};

%% Plot Resistance
p1.MarkerSize=0.5;
p1.LineWidth=1.5;
res=(Sim.Data.Rmat{IndexTime});
res=res(threshold,threshold);
[j,i,~]=find(tril(Adj2));
wd=zeros(1,length(j));
for k=1:length(j)
    wd(k)=res(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
[j,i,~]=find(tril(Adj2));
wd2=zeros(1,length(j));

for k=1:length(j)
    wd2(k)=Graph.networkThreshold(i(k),j(k));
end
wd3=wd(logical(wd2));

clim=[min([Sim.Settings.Roff Sim.Settings.Ron]) max([Sim.Settings.Roff Sim.Settings.Ron])];
p1.EdgeCData=wd3;
colormap(currAx,flipud(gcurrmap));%flipud(gray);
colorbar(currAx);
caxis(currAx,clim);

%Highlight Electrodes:
highlight(p1,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
labelnode(p1,highlightElec,[new_electrodes(:).Name]); %need to make this automated.

title(['Resistance Graph View | T= ' num2str(IndexTime)]);

%Voltage at each Node: %17/05/19
f5=figure;
currAx=gca;
p2=plot(currAx,G);
% set(currAx,'Color',[0.35 0.35 0.35]); %change background
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p2.NodeColor='red';
p2.EdgeColor='black';
p2.NodeLabel={};

%%Plot Voltage (log10)
if ~drain_exist %if no drains
    vlist=log10(Sim.Data.Voltages{IndexTime});
else
    vlist=Sim.Data.Voltages{IndexTime};
end
p2.NodeCData=full(vlist(threshold));
p2.MarkerSize=3;
colormap(currAx,hot);
colorbar(currAx);
caxis([Sim.SimInfo.MinV Sim.SimInfo.MaxV]);

labelnode(p2,highlightElec,[new_electrodes(:).Name]);
title(['Voltage Graph View |T= ' num2str(IndexTime)]);

Explore.GraphView.currents=Sim.Data.Currents{IndexTime};
Explore.GraphView.resistance=Sim.Data.Rmat{IndexTime};
Explore.GraphView.Nodes=G.Nodes;
Explore.GraphView.Edges=G.Edges;
Explore.GraphView.ElectrodePosition=Sim.Electrodes.PosIndex;
Explore.GraphView.AdjMat=Adj;
Explore.GraphView.Graph=G;
end

