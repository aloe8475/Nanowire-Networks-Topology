%% Network View:
% This function plots the network view of current, voltage and resistance
% for the chosen Sim at the given timestamp (IndexTime)

function [f1, f2, Adj, NumEl, Explore]= network_view(Sim, IndexTime, NodeList)

Layout=Sim.SelLayout;

% These are all aspects of the network that are used to graph it.
x1=diag(Layout.X1);
x2=diag(Layout.X2);
y1=diag(Layout.Y1);
y2=diag(Layout.Y2);
X=full([x1' ; x2']); % X = Wires 'x' value
Y=full([y1' ; y2']); % Y = Wires 'y' value
[~,~,Cx]=find(Layout.CX); %CX = Junctions 'x' value
[~,~,Cy]=find(Layout.CY); % CY = Junctions 'y' value
Adj=triu(Layout.AdjMat); % Adjacency matrix
NumEl=height(Sim.Electrodes); %Number of electrodes

%Plot Network:
f1=figure;
currAx=gca; %current axis
plot(currAx,X,Y,'b'); % plot wires
hold on
scatter(currAx,Cx,Cy,2,'r'); %scatterplot junctions
hold(currAx,'on');

%PlotElectrodes
IdxEl=Sim.Electrodes.PosIndex(NodeList.Value(NodeList.Value<=NumEl)); %Find Electrode Junction Position
Xe=X(:,IdxEl);Ye=Y(:,IdxEl); %Find X and Y (wires) for electrode
Cxe=(x2(IdxEl)+x1(IdxEl))./2;Cye=(y1(IdxEl)+y2(IdxEl))./2; %Find X and Y (junctions) for electrode
text(currAx,Cxe-1.7,Cye+0.7,'Electrode'); %Write 'Electrode' where it is placed

%Plot Currents:
currs=triu(Sim.Data.Currents{IndexTime}); % Find values of currents
Imat=full(abs(currs)); %full current matrix (instead of sparse double) + absolute value
Cx=Layout.CX(Adj~=0); %'x' coordinates for junctions that are 1 in the Adj matrix
Cy=Layout.CY(Adj~=0);%'y' coordinates for junctions that are 1 in the Adj matrix
Ilist=Imat(Adj~=0); % current list for junctions that are 1 in the Adj matrix
I=linspace(0,Sim.SimInfo.MaxI,10*length(Ilist)); %generates 10*length(Ilist)) points between 0 and max current
cmap=jet(10*length(Ilist)); %creates jet colormap with 10*length(Ilist) colorvalues
c=interp1(I,cmap,full(Ilist)); % %creates interpolation table with colormap - not sure what this is
PlotNetworkAux(currAx,X,Y,Cx,Cy,'curr',c);
labels=strsplit(num2str(1:length(Ilist))); %all junctions
% text(Cx,Cy,labels,'HorizontalAlignment','left');%label each junction with its number
clim=[Sim.SimInfo.MinI Sim.SimInfo.MaxI]; %minimum and maximum currents

%colorbar
colormap(currAx,cmap);
colorbar(currAx);
caxis(currAx,clim);
title(['Current Flow Network View Timestamp ' num2str(IndexTime)]);

%Plot Network for Figure 2
f2=figure;
currAx=gca; %current axis
plot(currAx,X,Y,'b'); % plot wires
hold on
scatter(currAx,Cx,Cy,2,'r'); %scatterplot junctions
hold(currAx,'on');

%PlotElectrodes
IdxEl=Sim.Electrodes.PosIndex(NodeList.Value(NodeList.Value<=NumEl)); %Find Electrode Junction Position
Xe=X(:,IdxEl);Ye=Y(:,IdxEl); %Find X and Y (wires) for electrode
Cxe=(x2(IdxEl)+x1(IdxEl))./2;Cye=(y1(IdxEl)+y2(IdxEl))./2; %Find X and Y (junctions) for electrode
text(currAx,Cxe-1.7,Cye+0.7,'Electrode'); %Write 'Electrode' where it is placed


%Plot Resistance
Rmat=triu(Sim.Data.Rmat{IndexTime});
Cx=Layout.CX(Adj~=0);
Cy=Layout.CY(Adj~=0);
Rlist=Rmat(Adj~=0);
R=linspace(min([Sim.Settings.Roff Sim.Settings.Ron]),max([Sim.Settings.Roff Sim.Settings.Ron]),10*length(Rlist));
cmap=flipud(copper(10*length(Rlist)));
c=interp1(R,cmap,full(Rlist));
PlotNetworkAux(currAx,X,Y,Cx,Cy,'curr',c);
clim=[min([Sim.Settings.Roff Sim.Settings.Ron]) max([Sim.Settings.Roff Sim.Settings.Ron])];
%colorbar
colormap(currAx,cmap);
colorbar(currAx);
caxis(currAx,clim);

title(['Resistance Network View Timestamp ' num2str(IndexTime)]);

%Save struct
Explore.NetworkView.currents=Ilist; %save currents at each junction at the IndexTime
Explore.NetworkView.resistance=Rlist; %save currents at each junction at the IndexTime
Explore.NetworkView.junctions=labels;
Explore.NetworkView.ElectrodePosition=IdxEl;

    end 