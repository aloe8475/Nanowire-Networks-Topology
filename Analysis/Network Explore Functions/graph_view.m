%% Graph View
% This function plots graph parameters such as current, voltage and
% resistance for the chosen Sim at the given timestamp (IndexTime)
function [f3, f4, f5, G, Adj, Explore,  highlightElec, new_electrodes]= graph_view(Sim,IndexTime,Explore,G, threshold_network)
if threshold_network=='t'
    return
    fprintf('Error in graph_view - you should not be seeing this');
end 
%Plot Graph
Adj=Sim.Data.AdjMat{IndexTime};
f3=figure;
currAx=gca;
p=plot(currAx,G);
% set(currAx,'Color',[0.35 0.35 0.35]); %change background
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p.NodeColor='red';
p.EdgeColor='white';
p.NodeLabel={};

%% Plot Currents
p.MarkerSize=1.5;
p.LineWidth=1.5;
currs=(abs(Sim.Data.Currents{IndexTime}));
[j,i,~]=find(tril(Adj));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end
clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
p.EdgeCData=cc;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);

%Highlight Electrodes:
for i=1:size(Sim.Electrodes.PosIndex,1)
    new_electrodes(i).PosIndex=Sim.Electrodes.PosIndex(i);
    new_electrodes(i).Name=Sim.Electrodes.Name(i);
end

highlightElec={new_electrodes.PosIndex};
highlightElec=cell2num(highlightElec);
highlight(p,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
labelnode(p,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
title(['Current Graph View Timestamp ' num2str(IndexTime)]);


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
[j,i,~]=find(tril(Adj));
wd=zeros(1,length(j));
for k=1:length(j)
    wd(k)=res(i(k),j(k));
end
clim=[min([Sim.Settings.Roff Sim.Settings.Ron]) max([Sim.Settings.Roff Sim.Settings.Ron])];
p1.EdgeCData=wd;
colormap(currAx,flipud(gcurrmap));%flipud(gray);
colorbar(currAx);
caxis(currAx,clim);

%Highlight Electrodes:
highlight(p1,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
labelnode(p1,highlightElec,[new_electrodes(:).Name]); %need to make this automated.

title(['Resistance Graph View Timestamp ' num2str(IndexTime)]);

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
vlist=Sim.Data.Voltages{IndexTime};
p2.NodeCData=full(vlist);
p2.MarkerSize=3;
colormap(currAx,hot);
colorbar(currAx);
caxis([Sim.SimInfo.MinV Sim.SimInfo.MaxV]);

labelnode(p2,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
title(['Voltage Graph View Timestamp ' num2str(IndexTime)]);

Explore.GraphView.currents=Sim.Data.Currents{IndexTime};
Explore.GraphView.resistance=Sim.Data.Rmat{IndexTime};
Explore.GraphView.Nodes=G.Nodes;
Explore.GraphView.Edges=G.Edges;
Explore.GraphView.ElectrodePosition=Sim.Electrodes.PosIndex;
end 