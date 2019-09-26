%% Graph Theory Explore
% This function plots graph theory parameters overlayed on graph view of
% currents for the chosen Sim at the given timestamp (IndexTime)

function [f6, f7, f8, f9, f10, f11, f12,f13, Explore, sourceElec, drainElec]= graph_theory_explore_threshold(Sim,G,Adj, Adj2, IndexTime,threshold,threshold_network, Explore, Graph, highlightElec, new_electrodes,node_indices,drain_exist,source_exist,network_load)
%% Find Source and Drain Electrodes:

%NOTE: This only works with 1 source and 1 drain:
if source_exist & ~drain_exist
    electrodes_cell(1)=new_electrodes(1).Name
elseif ~source_exist & drain_exist
    electrodes_cell(1)=new_electrodes(2).Name
else
    %this works as many sources/drains as required
    for i = 1:length(new_electrodes)
        electrodes_cell{i}=new_electrodes(i).Name;%% ALON TO FIX  06/08/19
    end
end
if ~isempty(new_electrodes)
    sourceIndex = find(contains(electrodes_cell,'Source')); %find index of electrodes that are source electrodes
    drainIndex = find(contains(electrodes_cell,'Drain')); %find index of electrodes that are drain electrodes
    
    if drain_exist
        drainElec=highlightElec(drainIndex);
    else
        drainElec=[];
    end
    
    if source_exist
        sourceElec=highlightElec(sourceIndex); %change to show path from different electrodes
    else
        sourceElec=[];
    end
    
    if ~isempty(sourceElec) & ~isempty(drainElec)
        noPath=0; %if we have both source and drain - we do have a path, so noPath = 0;
    elseif isempty(sourceElec) | isempty(drainElec)
        noPath=1;
    end
    
    for i =1:length(sourceElec)
        if source_exist
            source(i)=sourceElec(i); %choose first electrode if there are more than 1
        end
    end
    for i =1:length(drainElec)
        if drain_exist
            drain(i)=drainElec(i);%choose first electrode if there are more than 1
        end
    end
else
    noPath=1;
end
%% Participant Coefficients
if threshold_network~='t'
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
if network_load=='a'
    Adj=(Sim.Data.AdjMat{IndexTime});%we need to keep a copy of the original Adj matrix (unthresholded) to find all the currents - this is the NON BINARISED Adj matrix
else
    Adj=Sim.SelLayout.AdjMat;
end
Adj2=Adj(threshold,threshold);
currs=(abs(Sim.Data.Currents{IndexTime}));
currs=currs(threshold,threshold);
[j,i,~]=find(tril(Adj2));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
[j,i,~]=find(tril(Adj2));
cc2=zeros(1,length(j));

for k=1:length(j)
    cc2(k)=Graph.networkThreshold(i(k),j(k));
end

%Find CURRENTS in network
cc3=cc(logical(cc2));
if network_load=='a'
    clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
else
    clim=[min(min(Sim.Data.Currents{IndexTime})) max(max(Sim.Data.Currents{IndexTime}))];
end
p3.EdgeCData=cc3;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);
hold on

% Plot Participant Coefficients
p3.NodeCData=Graph.Ci(threshold);
p_ranks=Graph.P(threshold);
if ~noPath
    edges2 = linspace(min(p_ranks),max(p_ranks),7);
    bins2 = discretize(p_ranks,edges2);
    p3.MarkerSize=bins2;
    
    %Label Source
    
    
    labelnode(p3,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
end
%Need to figure out how to change colormap for Nodes seperately.
text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P(threshold))) ' (large dot)']);
title(['Participant Coeff, Thresholded | T= ' num2str(IndexTime)]);

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
if network_load=='a'
    Adj=(Sim.Data.AdjMat{IndexTime});%we need to keep a copy of the original Adj matrix (unthresholded) to find all the currents - this is the NON BINARISED Adj matrix
else
    Adj=Sim.SelLayout.AdjMat;
end
Adj2=Adj(threshold,threshold);
currs=(abs(Sim.Data.Currents{IndexTime}));
currs=currs(threshold,threshold);
[j,i,~]=find(tril(Adj2));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
[j,i,~]=find(tril(Adj2));
cc2=zeros(1,length(j));

for k=1:length(j)
    cc2(k)=Graph.networkThreshold(i(k),j(k));
end

%Find currents in network
cc3=cc(logical(cc2));
if network_load=='a'
    clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
else
    clim=[min(min(Sim.Data.Currents{IndexTime})) max(max(Sim.Data.Currents{IndexTime}))];
end
p4.EdgeCData=cc3;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);

hold on

p4.NodeCData=Graph.Ci(threshold);
mod_ranks=Graph.MZ(threshold);
if ~noPath
    edges3 = linspace(min(mod_ranks),max(mod_ranks),7);
    bins3 = discretize(mod_ranks,edges3);
    p4.MarkerSize=bins3;
    
    %Label Source
    labelnode(p4,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
end

text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ(threshold))) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);
title(['Within Module z-Score, Thresholded | T= ' num2str(IndexTime)]);

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
if network_load=='a'
    Adj=(Sim.Data.AdjMat{IndexTime});%we need to keep a copy of the original Adj matrix (unthresholded) to find all the currents - this is the NON BINARISED Adj matrix
else
    Adj=Sim.SelLayout.AdjMat;
end
Adj2=Adj(threshold,threshold);
currs=(abs(Sim.Data.Currents{IndexTime}));
currs=currs(threshold,threshold);
[j,i,~]=find(tril(Adj2));
cc=zeros(1,length(j));
for k=1:length(j)
    cc(k)=currs(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
[j,i,~]=find(tril(Adj2));
cc2=zeros(1,length(j));

for k=1:length(j)
    cc2(k)=Graph.networkThreshold(i(k),j(k));
end

%Find currents in network
cc3=cc(logical(cc2));
if network_load=='a'
    clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
else
    clim=[min(min(Sim.Data.Currents{IndexTime})) max(max(Sim.Data.Currents{IndexTime}))];
end
p5.EdgeCData=cc3;
colormap(currAx,gcurrmap);%gcurrmap
colorbar(currAx);
caxis(currAx,clim);

hold on
if ~noPath
    Graph.DEG_threshold=Graph.DEG(threshold);
    edges = linspace(min(Graph.DEG_threshold),max(Graph.DEG_threshold),7);
    bins = discretize(Graph.DEG_threshold,edges);
    p5.MarkerSize = bins;
    p5.NodeColor='r';
    title(['Degree Size, Thresholded | T= ' num2str(IndexTime)]);
    
    highlight(p5,highlightElec,'NodeColor','green'); %change simulation number
    labelnode(p5,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
end
text(-5,-6.2,'Min Degrees - 1 (small dot) | Max Degrees - 44 (large dot)');

%% Shortest Path (Distance)
d = distances(G); %calculate all the shortest path distance across all node pairs in the network.

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

%Plot all paths from current Electrode

f10=figure;
if ~noPath
    currAx=gca;
    p6=plot(currAx,G);
    p6.NodeLabel=d(source(1),:); %label the shortest path from the source;
    p6.NodeCData=d(source(1),:);
    highlight(p6,source(1),'MarkerSize',8);
    colormap(currAx,jet);%gcurrmap
    colorbar(currAx);
    title(['Path Distances from Source Electrode, Thresholded | T=' num2str(IndexTime)]);
    %% Show shortest path from source to drain: %27/05/19
    if drain_exist & source_exist
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
        if network_load=='a'
            Adj=(Sim.Data.AdjMat{IndexTime});%we need to keep a copy of the original Adj matrix (unthresholded) to find all the currents - this is the NON BINARISED Adj matrix
        else
            Adj=Sim.SelLayout.AdjMat;
        end
        Adj2=Adj(threshold,threshold);
        currs=(abs(Sim.Data.Currents{IndexTime}));
        currs=currs(threshold,threshold);
        [j,i,~]=find(tril(Adj2));
        cc=zeros(1,length(j));
        for k=1:length(j)
            cc(k)=currs(i(k),j(k));
        end
        
        % extract lower triangular part of Adjacency matrix of network
        [j,i,~]=find(tril(Adj2));
        cc2=zeros(1,length(j));
        
        %Find current in thresholded network:
        for k=1:length(j)
            cc2(k)=Graph.networkThreshold(i(k),j(k));
        end
        
        
        cc3=cc(logical(cc2));
        if network_load=='a'
            clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI];
        else
            clim=[min(min(Sim.Data.Currents{IndexTime})) max(max(Sim.Data.Currents{IndexTime}))];
        end
        p8.EdgeCData=cc3;
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
        if ~noPath
            p7.NodeLabel={};
            highlight(p7,highlightElec,'NodeColor','green','Marker','o'); %change simulation number
            [dist,path,pred]=graphshortestpath(Adj2,source(1),drain(1),'Directed','false');
            if length(sourceElec)>1
                [dist2,path2,pred2]=graphshortestpath(Adj2,source2,drain2,'Directed','false');
                highlight(p7,path2,'EdgeColor','cyan','LineWidth',6,'LineStyle','-');
            end
            highlight(p7,path,'EdgeColor','cyan','LineWidth',6,'LineStyle','-');
        end
        title(['Shortest Path, Thresholded + Overlayed on Current | T=' num2str(IndexTime)]);
    end
else
    f11=figure;
end


%% Communicability

f12=figure;
currAx=gca;
p9=plot(currAx,G);

Adj=Adj(threshold,threshold);
Graph.COMM=Graph.COMM(threshold,threshold);
% extract lower triangular part of Adjacency matrix of Graph.COMM
[j,i,~]=find(tril(Adj));
com=zeros(1,length(j));

for k=1:length(j)
    com(k)=Graph.COMM(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
if network_load=='a'
    Adj2=(Sim.Data.AdjMat{IndexTime});%we need to keep a copy of the original Adj matrix (unthresholded) to find all the currents - this is the NON BINARISED Adj matrix
else
    Adj2=Sim.SelLayout.AdjMat;
end
Adj2=Adj2(threshold,threshold);
[j,i,~]=find(tril(Adj2));
com2=zeros(1,length(j));

for k=1:length(j)
    com2(k)=Graph.networkThreshold(i(k),j(k));
end

%Find Graph.COMM in network
com3=com(logical(com2));

p9.EdgeCData=com3;%log10(com3);
p9.MarkerSize=2;
colormap jet
colorbar
% labelnode(p9,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
if ~noPath
    labelnode(p9,highlightElec,[new_electrodes(:).Name]);
    
    highlight(p9,highlightElec,'NodeColor','green','Marker','o'); %change simulation number
end
%Overlay Shortest Path
if ~noPath
    hold on
    p10=plot(currAx,G);
    set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
    p10.NodeColor='green';
    p10.Marker='none';
    p10.EdgeColor='green';
    p10.LineStyle='none';
    p10.NodeLabel={};
    highlight(p10,highlightElec,'NodeColor','green','Marker','o'); %change simulation number
    [dist,path,pred]=graphshortestpath(Adj,source(1),drain(1),'Directed','false');
    if length(sourceElec)>1
        [dist2,path2,pred2]=graphshortestpath(Adj,source(2),drain(2),'Directed','false');
        highlight(p10,path2,'EdgeColor','[0.95, 0.95, 0.95]','LineWidth',6,'LineStyle','-');
    end
    highlight(p10,path,'EdgeColor','[0.95, 0.95, 0.95]','LineWidth',6,'LineStyle','-');
    title(['Communicability, Thresholded + Overlayed w Shortest Path | T=' num2str(IndexTime)]);
end
%% Clustering Coefficient

%Cluster Analysis:
f13=figure;
currAx=gca;
p11=plot(currAx,G);
p11.MarkerSize = 4;
% set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
p11.NodeCData=Graph.Ci(threshold);
% Ci_ranks=Graph.Ci(threshold);
% edges3 = linspace(min(Ci_ranks),max(Ci_ranks),7);
% bins3 = discretize(Ci_ranks,edges3);
% p11.MarkerSize=bins3;

colormap hsv(6)
% colorbar
if ~noPath
    labelnode(p11,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
    labelnode(p11,highlightElec,[new_electrodes(:).Name]);
end
if ~noPath
    %Overlay Shortest Path
    hold on
    p12=plot(currAx,G);
    set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
    p12.NodeColor='green';
    p12.Marker='none';
    p12.EdgeColor='green';
    p12.LineStyle='none';
    p12.NodeLabel={};
    % highlight(p12,highlightElec,'NodeColor','green','Marker','o'); %change simulation number
    [dist,path,pred]=graphshortestpath(Adj,source(1),drain(1),'Directed','false');
    if length(sourceElec)>1
        [dist2,path2,pred2]=graphshortestpath(Adj,source(2),drain(2),'Directed','false');
        highlight(p12,path2,'EdgeColor','[0.95, 0.95, 0.95]','LineWidth',6,'LineStyle','-');
    end
    highlight(p12,path,'EdgeColor','[0.95, 0.95, 0.95]','LineWidth',6,'LineStyle','-');
    title(['Clustering, Thresholded + Overlayed w Shortest Path | T=' num2str(IndexTime)]);
end

Explore.GraphView.Distances.Values=d;
Explore.GraphView.Distances.Avg=avgD;
Explore.GraphView.Distances.Std=stdD;
Explore.GraphView.Distances.Median=medianD;
if ~noPath
    Explore.GraphView.Distances.DistancesFromSource=d(sourceElec,:);
    if source_exist
        Explore.Electrodes.Source=sourceElec;
    else
        Explore.Electrodes.Source=[];
    end
    if drain_exist
        Explore.Electrodes.Drain=drainElec;
    else
        Explore.Electrodes.Drain=[]
    end
    if drain_exist & source_exist %if both source and drain - we have shortest paths
        Explore.GraphView.Distances.ShortestPath=path;
    end
else
    sourceElec=[];
    drainElec=[];
end
end