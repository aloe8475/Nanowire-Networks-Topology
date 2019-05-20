%---------------------------
% 22/04/19
% Network Explore Code v1.1
%
% Author: Alon Loeffler
%
% Changelog:
% 13/05/19 - Added support for simulation files with only 1 simulation
% --------------------------
dbstop if error
%% Load Data
load_data_question=lower(input('Load network data, Analysis Data Only or None? N - None, D - Network Data, A - Analysis Data\n','s'));

if load_data_question=='d'
    clearvars -except load_data_question
    close all
    %load network data
    [network, network_load, simulations,sim_loaded,numNetworks, explore_network]= load_data();
elseif load_data_question=='a'
    clear LDA_Analysis
    close all
    %Load previous LDA analysis data
    networkNum=input(['Which Network # do you want to load? 1 - ' num2str(length(network)) '\n']);
    simNum=input(['Which Simulation # do you want to load? 1 - '  num2str(length(network(networkNum).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    [LDA_Analysis(simNum)] = load_LDA_data();
end

%% Choose Simulations
%-------
%Choose Network and Simulation for training or exploring
if load_data_question~='a'
    if explore_network=='t'
        networkNum=input(['Which Network # do you want to select for Training? 1 - ' num2str(length(network)) '\n']);
        simNum=input(['Which Simulation # do you want to select for Training? 1 - '  num2str(length(network(networkNum).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    else
        networkNum=input(['Which Network # do you want to explore? 1 - ' num2str(length(network)) '\n']);
        simNum=input(['Which Simulation # do you want to explore? 1 - '  num2str(length(network(networkNum).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    end
end

if size(simulations,1)==1
    currentSim=simulations{simNum};
else
    currentSim=simulations(simNum);
end

fprintf(['Simulation: ' network(networkNum).Name currentSim.Name ' selected \n\n']);

%% Analysis:
i = 1;
while i == 1
    %Choose Analysis to perform
    if explore_network=='t'
        analysis_type=lower(input('Which analysis would you like to perform? G - graph, E - Explore Network,L - LDA, N - none \n','s'));
    elseif explore_network=='e' % we don't want to allow LDA if just exploring
        analysis_type=lower(input('Which analysis would you like to perform? G - graph, E - Explore Network, N - none \n','s'));
    end
    if analysis_type=='g'
    %% Graph Analysis
        [Graph, threshold_network]=graph_analysis(network(networkNum),network_load,currentSim,[]);
        %% Save Graph analysis data
        i=i+1;
    elseif analysis_type=='e'
        %% Exploratary analysis of simulation
        Explore=explore_simulation(currentSim,network,network_load,simNum);
        %% Saving Explore
        save_state=lower(input('Would you like to save the Exploration Analysis? y or n \n','s'));
        if save_state=='y'
            save_explore(Explore,network(networkNum),network_load);
            i=i+1; %get out of while loop this loop finishes
        end
        i=i+1;
    elseif analysis_type=='l'
        %% LDA Analysis
        LDA_Analysis=lda_analysis(currentSim,network,network_load,simNum);
        i=i+1;
    end
    if analysis_type=='l' || analysis_type=='n'
        %% Plotting LDA
        plot_state=lower(input('Would you like to plot LDA Analysis? y or n \n','s'));
        if plot_state=='y'
            plot_LDA(LDA_Analysis(simNum),simNum,network(networkNum).Name);
            i=i+1;
        end
        
        %% Apply LDA Training to testing data (different simulation)
        if analysis_type~='n'
            apply_LDA=lower(input('Would you like to apply the loaded LDA analysis to another simulation? y or n \n','s'));
            if apply_LDA=='y'
                [LDA_Analysis, simulationChoice]=lda_apply_func(numNetworks,network,LDA_Analysis,simNum);
            end
            i=i+1; %get out of while loop when this loop finishes
        elseif analysis_type~='l' &&  analysis_type~='g' && analysis_type~='n' && analysis_type~='e'
            fprintf('Please type either G, L or N only \n');
        end
        
        if analysis_type=='g' || analysis_type=='n'
            %% Plotting Graph
            plot_state=lower(input('Would you like to plot Graph Analysis? y or n \n','s'));
            if plot_state=='y'
                Graph=plot_graph(Graph,network(networkNum),network_load,currentSim,sim_loaded);
            end
            i=i+1;
            %% Saving Graph
            save_state=lower(input('Would you like to save the Graph Analysis? y or n \n','s'));
            if save_state=='y'
                save_graph(Graph,network(networkNum),network_load);
            end
            i=i+1;
        end
        %% Insert Further Analysis Below
        % -------------------------------
        % -------------------------------
    end
end

%--------------------------------------------------------------------------


%% FUNCTIONS

function [network, network_load, simulations, sim_loaded, numNetworks, explore_network] = load_data()

%% Load Data
%Ask to load Zdenka or Adrian:
network_load=lower(input('Which Network do you want to analyse? Z - Zdenka, A - Adrian \n','s'));

if strcmp(network_load,'a')
    %Get current network - Adrian
    [network,sim_loaded, explore_network, numNetworks]=Load_Adrian_Code();
    %unpack simulation data into simulation variable
    if sim_loaded==1
        if explore_network=='t' %if we have training and testing simulations
            tempSim=network.Simulations{2};
            tempSim=num2cell(tempSim);
            network.Simulations(2) = [];
            network.Simulations=[network.Simulations tempSim];
            fprintf(['Your Training Simulation is Simulation 1 \n']);
            fprintf(['Your Testing Simulations are Simulations 2 - ' num2str(length(network.Simulations)) '\n']);
            
            fprintf('\n -------------------------- \nStart Analysis: \n');
            
        end
        for i = 1:length(network.Simulations)
            simulations(i)=network.Simulations(i);
        end
    else
        simulations=network.Simulations;
    end
    cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis');
elseif strcmp(network_load,'z')
    %Get network - Zdenka:
    % D:\alon_\Research\POSTGRAD\PhD\CODE\Zdenka's Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data
    network=Load_Zdenka_Code();
    cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis');
end
end

function Explore = explore_simulation(Sim,network,network_load,simNum)
[NodeList.String,NodeList.UserData]=GetNodeList(Sim);
NodeList.Value=1;

%% Network View:
Layout=Sim.SelLayout;

% These are all aspects of the network that are used to graph it.
x1=diag(Layout.X1);
x2=diag(Layout.X2);
y1=diag(Layout.Y1);
y2=diag(Layout.Y2);
X=full([x1' ; x2']); % X = Wires 'x' value
Y=full([y1' ; y2']); % Y = Wires 'y' value
Gr=Layout.SelGraph; % graph
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

%Choose a time to Explore Simulation:
IndexTime=input(['What Timestamp do you want to analyse? 1-' num2str(size(Sim.Data,1)) '\n']); %CHOOSE TIMESTAMP

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
cmap=flipud(gray(10*length(Rlist)));
c=interp1(R,cmap,full(Rlist));
PlotNetworkAux(currAx,X,Y,Cx,Cy,'curr',c);
clim=[min([Sim.Settings.Roff Sim.Settings.Ron]) max([Sim.Settings.Roff Sim.Settings.Ron])];
%colorbar
colormap(currAx,cmap);
colorbar(currAx);
caxis(currAx,clim);

title(['Resistance Network View Timestamp ' num2str(IndexTime)]);


%% Graph View

%Plot Graph
Layout=Sim.SelLayout;
G=Layout.SelGraph;
Adj=(Layout.AdjMat);
f3=figure;
currAx=gca;
p=plot(currAx,G);
set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p.NodeColor='red';
p.EdgeColor='white';
p.NodeLabel={};

%Plot Currents (log10)
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
set(currAx,'Color',[0.35 0.35 0.35]);
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p1.NodeColor='red';
p1.EdgeColor='white';
p1.NodeLabel={};

%Plot Resistance
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
set(currAx,'Color',[0.35 0.35 0.35]);
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p2.NodeColor='red';
p2.EdgeColor='white';
p2.NodeLabel={};

%Plot Voltage (log10)
vlist=Sim.Data.Voltages{IndexTime};
p2.NodeCData=full(vlist);
p2.MarkerSize=3;
colormap(currAx,hot);
colorbar(currAx);
caxis([Sim.SimInfo.MinV Sim.SimInfo.MaxV]);

labelnode(p2,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
title(['Voltage Graph View Timestamp ' num2str(IndexTime)]);

%% Overlay Graph Theory:
[Graph, threshold_network]=graph_analysis(network,network_load,Sim,IndexTime);

%% Participant Coefficients

f6=figure;
currAx=gca;
p3=plot(currAx,G);
set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p3.NodeColor='red';
p3.EdgeColor='white';
p3.NodeLabel={};

%Plot Currents
p3.MarkerSize=1.5;
p3.LineWidth=1.5;
currs=(abs(Sim.Data.Currents{IndexTime}));
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
if threshold_network=='y'
    p3.NodeCData=Graph.Ci(threshold);
    p_ranks=Graph.P(threshold);
else
    p3.NodeCData=Graph.Ci;
    p_ranks=Graph.P;
end
edges2 = linspace(min(p_ranks),max(p_ranks),7);
bins2 = discretize(p_ranks,edges2);
p3.MarkerSize=bins2;

%Label Source
labelnode(p3,highlightElec,[new_electrodes(:).Name]); %need to make this automated.

%Need to figure out how to change colormap for Nodes seperately.

if threshold_network=='y'
    text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P(threshold))) ' (large dot)']);
else
    text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P)) ' (large dot)']);
end
title(['Participant Coefficient Analysis Timestamp ' num2str(IndexTime)]);

%% Modular z-Score:
f7=figure;
currAx=gca;
p4=plot(currAx,G);
set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p4.NodeColor='red';
p4.EdgeColor='white';
p4.NodeLabel={};

%Plot Currents
p4.MarkerSize=1.5;
p4.LineWidth=1.5;
currs=(abs(Sim.Data.Currents{IndexTime}));
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

if threshold_network=='y'
    p4.NodeCData=Graph.Ci(threshold);
    mod_ranks=Graph.MZ(threshold);
else
    p4.NodeCData=Graph.Ci;
    mod_ranks=Graph.MZ;
end
edges3 = linspace(min(mod_ranks),max(mod_ranks),7);
bins3 = discretize(mod_ranks,edges3);
p4.MarkerSize=bins3;

%Label Source
labelnode(p4,highlightElec,[new_electrodes(:).Name]); %need to make this automated.

if threshold_network=='y'
    text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ(threshold))) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);
else
    text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ)) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);
end
title(['Within Module Degree z-Score Analysis Timestamp ' num2str(IndexTime)]);

%% Connectivity
f8=figure;
currAx=gca;
p5=plot(currAx,G);
set(currAx,'Color',[0.35 0.35 0.35]);% change background color to gray
set(gcf, 'InvertHardCopy', 'off'); %make sure to keep background color
p5.NodeColor='red';
p5.EdgeColor='white';
p5.NodeLabel={};

%Plot Currents
p5.MarkerSize=1.5;
p5.LineWidth=1.5;
currs=(abs(Sim.Data.Currents{IndexTime}));
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

if threshold_network=='y'
    Graph.DEG_threshold=Graph.DEG(threshold);
    edges = linspace(min(Graph.DEG_threshold),max(Graph.DEG_threshold),7);
    bins = discretize(Graph.DEG_threshold,edges);
else
    edges = linspace(min(Graph.DEG),max(Graph.DEG),7);
    bins = discretize(Graph.DEG,edges);
end

p5.MarkerSize = bins;
p5.NodeColor='r';
title(['Degree Size Timestamp ' num2str(IndexTime)]);
highlight(p5,highlightElec,'NodeColor','green'); %change simulation number
labelnode(p5,highlightElec,[new_electrodes(:).Name]); %need to make this automated.
text(-5,-6.2,'Min Degrees - 1 (small dot) | Max Degrees - 44 (large dot)');

%Shortest Path (Distance)
d = distances(G); %calculate all the shortest path distance across all node pairs in the network.

%% Save
%Save Variables
Explore.IndexTime=IndexTime;
Explore.Name=Sim.Name;
Explore.NetworkView.currents=Ilist; %save currents at each junction at the IndexTime
Explore.NetworkView.resistance=Rlist; %save currents at each junction at the IndexTime
Explore.NetworkView.junctions=labels;
Explore.NetworkView.ElectrodePosition=IdxEl;
Explore.GraphView.currents=Sim.Data.Currents{IndexTime};
Explore.GraphView.resistance=Sim.Data.Rmat{IndexTime};
Explore.GraphView.Nodes=G.Nodes;
Explore.GraphView.Edges=G.Edges;
Explore.GraphView.ElectrodePosition=Sim.Electrodes.PosIndex;
Explore.GraphView.Distances=d;

%% Save Plots
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\Explore Analysis\';
network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems

if threshold_network~='y'
    saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'eps');
    saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'eps');
    saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Currents_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Currents_Timestamp' num2str(IndexTime)],'eps');
    saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Resistance_Timestamp' num2str(IndexTime)],'eps');
    saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Voltage_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Voltage_Timestamp' num2str(IndexTime)],'eps');
    saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'eps');
    saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'eps');
    saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'eps');
elseif threshold_network=='y'
    saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'eps');
    saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'eps');
    saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Currents_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Currents_Timestamp' num2str(IndexTime)],'eps');
    saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Resistance_Timestamp' num2str(IndexTime)],'eps');
    saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Voltage_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Voltage_Timestamp' num2str(IndexTime)],'eps');
    saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'eps');
    saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'eps');
    saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'jpg');
    saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'eps');
    
end
end

function save_explore(Explore,network,network_load)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Network Exploration Analysis\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_Exploration_Analysis_' date],'Explore');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'Exploration_Analysis' num2str(Explore.IndexTime) '_Timestamp_' date],'Explore');
end
end

function LDA_Analysis=lda_analysis(currentSim,network,network_load,simNum)
%% LDA Analysis
%THIS IS WHERE YOU CHANGE THE VARIABLES FOR LDA
drain1=[full(currentSim.Data.IDrain1)]; %full(simulations(2).Data.IDrain1)];
drain2=[full(currentSim.Data.IDrain2)]; %full(simulations(2).Data.IDrain2)];
%Current source
source1=[full(currentSim.Data.ISource1)]; %full(simulations(2).Data.ISource1)];
source2=[full(currentSim.Data.ISource2)]; %full(simulations(2).Data.ISource2)];
%Voltage source
source1Voltage=full(currentSim.Data.VSource1);
source2Voltage=full(currentSim.Data.VSource2);
Input = [source1 source2];
Output = [drain1 drain2]; %OUTPUT is column (variable) x row (observations)
Target = source1Voltage(:,1)>0.001;% uncomment for current: [source1(:,1)> 1.0e-04 *0.001 & source2<0];%TARGET is the classifier we expect

%visualise drain1 and drain2
figure
plot(drain1); hold on
plot(drain2);

% Calculate linear discriminant coefficients - finding a line that
% differentiate Output & Target
LDA_Analysis(simNum).W = LDA(Output,Target);

% % Calulcate linear scores for training data - what are the loadings for
% that line - how would you get to 'x' from that line
LDA_Analysis(simNum).L = [ones(size(Output,1),1) Output] * LDA_Analysis(simNum).W';
%
% % Calculate class probabilities
LDA_Analysis(simNum).P = exp(LDA_Analysis(simNum).L) ./ repmat(sum(exp(LDA_Analysis(simNum).L),2),[1 2]);

%Save all variables into struct
LDA_Analysis(simNum).drain1=drain1;
LDA_Analysis(simNum).drain2=drain2;
LDA_Analysis(simNum).source1=source1;
LDA_Analysis(simNum).source2=source2;
LDA_Analysis(simNum).Output=Output;
LDA_Analysis(simNum).Input=Input;
LDA_Analysis(simNum).Target=Target;

LDA_Analysis(simNum).normalisedOutput=LDA_normalise(Output);
LDA_Analysis(simNum).normalisedInput=LDA_normalise(Input);
LDA_Analysis(simNum).normalisedW=LDA(LDA_Analysis(simNum).normalisedOutput,Target);
LDA_Analysis(simNum).normalisedL=[ones(size(LDA_Analysis(simNum).normalisedOutput,1),1) LDA_Analysis(simNum).normalisedOutput] * LDA_Analysis(simNum).normalisedW';
LDA_Analysis(simNum).normalisedP=exp(LDA_Analysis(simNum).normalisedL) ./ repmat(sum(exp(LDA_Analysis(simNum).normalisedL),2),[1 2]);

LDA_Analysis(simNum).TypeOfData='Training';
clear drain1 drain2 source1 source2 Input Target Output
%% Saving LDA
save_state=lower(input('Would you like to save the LDA Analysis? y or n \n','s'));
if save_state=='y'
    save_LDA(LDA_Analysis(simNum),network(networkNum),network_load,sim_loaded);
end
end

function [LDA_Analysis, simulationChoice]=lda_apply_func(numNetworks,network,LDA_Analysis,simNum)
if numNetworks>1 %if we have two networks, offer to test second network
    networkNum2=input(['Which Network # do you want to select your Simulation from ? 1 - ' num2str(length(network)) '\n']);
else
    networkNum2=1;
end
if numNetworks>1 %if we have more than 1 simulation, they can input 1 as an option
    simulationChoice=input(['Which Simulation # do you want to apply LDA to? 1 - '  num2str(length(network(networkNum2).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
else % otherwise, can only input 2 or higher
    while 1
        simulationChoice=input(['Which Simulation # do you want to apply LDA to? 2 - '  num2str(length(network(networkNum2).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
        if simulationChoice >1 && simulationChoice <=length(network(networkNum2).Simulations)
            break;
        elseif simulationChoice ==1
            fprintf('Cannot Train and Test the same Simulation, please choose again \n');
        else
            fprintf(['Please Choose a number between 2 and ' num2str(length(network(networkNum2).Simulations))]);
        end
    end
    LDA_Analysis(simulationChoice).Output=[full(simulations{simulationChoice}.Data.IDrain1) full(simulations{simulationChoice}.Data.IDrain2)];
    LDA_Analysis(simulationChoice).Input=[full(simulations{simulationChoice}.Data.ISource1) full(simulations{simulationChoice}.Data.ISource2)];
    LDA_Analysis(simulationChoice).Target=[full(simulations{simulationChoice}.Data.VSource1 >0.001)];
    LDA_Analysis(simulationChoice).TypeOfData='Testing';
    [LDA_Analysis(simulationChoice).appliedP, LDA_Analysis(simulationChoice).appliedL, LDA_Analysis(simulationChoice).normalisedOutput, LDA_Analysis(simulationChoice).normalisedInput]=LDA_Apply(LDA_Analysis(simNum).normalisedW,LDA_Analysis(simulationChoice).Output, LDA_Analysis(simulationChoice).Input);
    plot_state2=lower(input('Would you like to plot the applied LDA results? y or n \n','s'));
    if plot_state2=='y'
        plot_LDA_Applied(LDA_Analysis(simulationChoice),simNum,simulationChoice,network(networkNum2).Name);
    end
end
end

function [LDA_Analysis] = load_LDA_data()
cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Data\LDA Analysis (Mac)');
waitfor(msgbox('Select the LDA Analysis saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the LDA Analysis saved data');
f=fullfile(PathName,FileName);
load(f);
cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis');
end

function plot_LDA(LDA_Analysis, simNum, networkName)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\LDA\LDA Training\';
%plot LDA
LDAf=figure;
sgtitle('LDA Classification Training');
subplot(4,1,1)
plot(LDA_Analysis.Input)
title('Source Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
plot(LDA_Analysis.Output)
title('Drain Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.P(:,2))
title('LDA class probability (P)');
xlabel('Timestamp (0.01sec)')
ylabel('Probability');
subplot(4,1,4)
plot(LDA_Analysis.Target)
title('LDA Target Classification (class labels) - Source1 > 0V')
xlabel('Timestamp(0.01sec)');
ylabel('Voltage > 0');

%plot Log LDA
LDAff=figure;
sgtitle('Log-Scale LDA Classification Training');
subplot(4,1,1)
semilogy(LDA_Analysis.Input) %log y axis plot
title('Source Electrodes');
ylabel('Log Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
semilogy(LDA_Analysis.Output) %log y axis plot
title('Drain Electrodes');
ylabel('Log Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.P(:,2))
title('LDA class probability (P)');
xlabel('Timestamp (0.01sec)')
ylabel('Probability');
subplot(4,1,4)
plot(LDA_Analysis.Target)
title('LDA Target Classification (class labels) - Source1 > 0V')
xlabel('Timestamp(0.01sec)');
ylabel('Voltage > 0');

LDAnormf=figure;
sgtitle('Normalised LDA Classification Training');
subplot(4,1,1)
plot(LDA_Analysis.normalisedInput)
title('Normalised Source Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
plot(LDA_Analysis.normalisedOutput)
title('Normalised Drain Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.normalisedP(:,2))
title('LDA class probability (P)');
xlabel('Timestamp (0.01sec)')
ylabel('Probability');
subplot(4,1,4)
plot(LDA_Analysis.Target)
title('LDA Target Classification (class labels) - Source1 > 0V')
xlabel('Timestamp(0.01sec)');
ylabel('Voltage > 0');

networkName(regexp(networkName,'[/:]'))=[]; %remove '/' character because it gives us saving problems

%Save
%Note; change the date of the simulation in the name if using a different
%set of simulations
saveas(LDAf,[save_directory num2str(networkName) 'Simulation' num2str(simNum) '_LDA_Classification_Training_Simulation' num2str(simNum)],'jpg'); % CHANGE DATE OF SIMULATION
saveas(LDAf,[save_directory num2str(networkName) 'Simulation' num2str(simNum) '_LDA_Classification_Training_Simulation' num2str(simNum)],'eps');
saveas(LDAff,[save_directory num2str(networkName) 'Simulation' num2str(simNum) '_logLDA_Classification_Training_Simulation' num2str(simNum)],'jpg');
saveas(LDAff,[save_directory num2str(networkName) 'Simulation' num2str(simNum) '_logLDA_Classification_Training_Simulation' num2str(simNum)],'eps');
saveas(LDAnormf,[save_directory num2str(networkName) 'Simulation' num2str(simNum) '_normalised_LDA_Classification_Training_Simulation' num2str(simNum)],'jpg');
saveas(LDAnormf,[save_directory num2str(networkName) 'Simulation' num2str(simNum) '_normalised_Classification_Training_Simulation' num2str(simNum)],'eps');

end

function plot_LDA_Applied(LDA_Analysis,simNum,simulationChoice,networkName)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\LDA\LDA Testing\';
appliedF=figure;
sgtitle('Electrode LDA Classification Testing');
subplot(4,1,1)
plot(LDA_Analysis.normalisedInput);
title('Normalised Source Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,2)
plot(LDA_Analysis.normalisedOutput);
title('Normalised Drain Electrodes');
ylabel('Current (A)')
xlabel('Timestamp (0.01sec)')
subplot(4,1,3)
plot(LDA_Analysis.appliedP(:,2));
title('Classification Probability');
ylabel('Probability')
xlabel('Timestamp (0.01sec)')
subplot(4,1,4)
plot(LDA_Analysis.Target(:,1));
title('LDA Target Classification (class labels) - Source1 > 0V');
ylabel('Voltage > 0')
xlabel('Timestamp (0.01sec)')

networkName(regexp(networkName,'[/:]'))=[]; %remove '/' character because it gives us saving problems

saveas(appliedF,[save_directory num2str(networkName) '_6MaySim_LDA_Classification_TrainingSim_' num2str(simNum) '_TestingSim' num2str(simulationChoice)],'jpg');
saveas(appliedF,[save_directory num2str(networkName) '_6MaySim_LDA_Classification_TrainingSim_' num2str(simNum) '_TestingSim' num2str(simulationChoice)],'eps');
end

function [Graph, threshold_network]=graph_analysis(network,network_load,currentSim,IndexTime)
%% Mac's Analysis: (Graph)
if strcmp(network_load,'z')%Zdenka Code:
    net_mat=network.adj_matrix; %symmetrical matrix
elseif strcmp(network_load,'a') %adrian code
    if isempty(IndexTime)
        IndexTime=input(['What Timestamp do you want to analyse? 1-' num2str(size(currentSim.Data,1)) '\n']); %CHOOSE TIMESTAMP
    end
    %this gives a resistance matrix for the network used for a chosen simulation at a specific timestamp
    threshold_network=input('Do you want to threshold the network using Resistance? \n','s');
    if threshold_network=='y'
        a= full(currentSim.Data.Rmat{IndexTime}); %Using Resistance
        a(a==5000)=1;
        a(a==5000000)=0;  %binarising the resistance
        net_mat=a;
    else
        net_mat=currentSim.SelLayout.AdjMat; %use standard adjacency matrix
    end
    %    net_mat=full(simulations(simNum).Data.AdjMat{IndexTime}); %this gives an adj matrix for the network used for a chosen simulation at a specific timestamp
end

%Global efficiency --> 1/characteristic path length, averaged over the whole network. An estimate of how integrated the network is.
Graph.GE = efficiency_bin(net_mat,0);

%Topological features will allow us to understand whether it was good or
%bad for classification.

%Local efficiency --> 1/characteristic path length, at each node
Graph.LE = efficiency_bin(net_mat,1);

%communicability --> an estimate of the ease with which each pair of nodes can connect in the network
Graph.COMM = expm(net_mat);

%modularity --> an estimate of how segregated the network is
[Graph.Ci,Graph.Q] = community_louvain(net_mat,1);
%The Ci term is the module assignment for each node
%The Q term is the 'quality' of the partition --> how modular the network is.
% -- this should tell us how well we can classify

%degree --> a count of how many edges are connected to each node
Graph.DEG = degrees_und(net_mat);

%participation coefficient --> an estimate of how integrative a node is
Graph.P = participation_coef(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'

%module degree z-score --> an estimate of how segregated a node is
Graph.MZ = module_degree_zscore(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'


%Network Density:
% This is defined, for a given set of nodes, as the number of actual edges
% divided by the number of potential edges.
% I.e., the percentage of possible connections that actually exist.


%save network matrix to graph struct
Graph.network=net_mat;
Graph.IndexTime=IndexTime;
end

function Graph=plot_graph(Graph, network,network_load, currentSim,sim_loaded)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Figures\Graph Analysis\';

IndexTime=Graph.IndexTime;
%visualise graph network:
threshold=Graph.DEG>1; %greater than 1 degree threshold
Graph.networkThreshold=Graph.network(threshold,threshold); %applying degree threshold
g=graph(Graph.networkThreshold);

f1=figure;
p1=plot(g);

%find electrodes:
node_indices=find(threshold==1); %find nodes with threshold == 1
for i=1:size(currentSim.Electrodes.PosIndex,1)
    if ~isempty(find(node_indices==currentSim.Electrodes.PosIndex(i)))
        new_electrodes(i).PosIndex=find(node_indices==currentSim.Electrodes.PosIndex(i));
        new_electrodes(i).Name=currentSim.Electrodes.Name(i);
    end
end

highlightElec={new_electrodes.PosIndex};
highlightElec=cell2num(highlightElec);
%highlight electrodes on graph
if sim_loaded==1
    highlight(p1,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
    labelnode(p1,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
end

f2=figure;
p2=plot(g);
%Closeness graph: (Unweighted)
ucc = centrality(g,'closeness');
p2.NodeCData=ucc;
colormap jet
colorbar
title(['Closeness Centrality Scores Timestamp ' num2str(IndexTime)])
if sim_loaded==1
    highlight(p2,highlightElec,'MarkerSize',7); %change simulation number
    labelnode(p2,[1:size(node_indices,2)],cellstr(num2str(node_indices'))); %label each node with original node number
    labelnode(p2,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
end

%Node size graph:
f3=figure;
p3=plot(g);
Graph.DEG_threshold=Graph.DEG(threshold);
edges = linspace(min(Graph.DEG_threshold),max(Graph.DEG_threshold),7);
bins = discretize(Graph.DEG_threshold,edges);
p3.MarkerSize = bins;
p3.NodeColor='r';
title(['Degree Size Timestamp ' num2str(IndexTime)]);
if sim_loaded==1
    highlight(p3,highlightElec,'NodeColor','green'); %change simulation number
    labelnode(p3,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
    labelnode(p3,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
end
text(-5,-6.2,'Min Degrees - 1 (small dot) | Max Degrees - 44 (large dot)');

%Both combined:
f4=figure;
p4=plot(g);
p4.MarkerSize = bins;
p4.NodeCData=ucc;
colormap jet
colorbar
labelnode(p4,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number

labelnode(p4,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
title(['Centrality Score and Degree Size Timestamp ' num2str(IndexTime)]);
text(-5.5,-6.2,'Min Degrees = 1 (small dot) | Max Degrees = 44 (large dot)');

%Histogram of degree distribution:
f5=figure;
h1=histogram(Graph.DEG_threshold);
title(['Distribution of Connectivity of Nodes Timestamp ' num2str(IndexTime)]);
xlabel('Number of Connections (Degrees)');
ylabel('Frequency');
Graph.avgDEG=mean(Graph.DEG(threshold));
Graph.stdDEG=std(Graph.DEG(threshold));
ylim([0 30]);
text(4,16,['Mean: ' num2str(Graph.avgDEG) ' | SD: ' num2str(Graph.stdDEG)]);

%Cluster Analysis:
f6=figure;
p6=plot(g);
p6.MarkerSize = 4;
p6.NodeCData=Graph.Ci(threshold);
labelnode(p6,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number

labelnode(p6,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
colormap hsv(6) %change number of colors here if there are more/less than 6 clusters
title(['Cluster Analysis ' num2str(IndexTime)]);

%Participant Coefficient Analysis:
f7=figure;
p7=plot(g);
p7.NodeCData=Graph.Ci(threshold);
p_ranks=Graph.P(threshold);
edges2 = linspace(min(p_ranks),max(p_ranks),7);
bins2 = discretize(p_ranks,edges2);
p7.MarkerSize=bins2;
colormap hsv(6)
labelnode(p7,[1:size(node_indices,2)],cellstr(num2str(node_indices'))); %label each node with original node number

labelnode(p7,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P(threshold))) ' (large dot)']);
title(['Participant Coefficient Analysis Timestamp ' num2str(IndexTime)]);

%Module Degree Z-Score:
f8=figure;
p8=plot(g);
p8.NodeCData=Graph.Ci(threshold);
mod_ranks=Graph.MZ(threshold);
edges3 = linspace(min(mod_ranks),max(mod_ranks),7);
bins3 = discretize(mod_ranks,edges3);
p8.MarkerSize=bins3;
colormap hsv(6)
labelnode(p8,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number

labelnode(p8,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable
text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ(threshold))) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);
title(['Module Degree z-Score Analysis Timestamp ' num2str(IndexTime)]);


% Communicability at different times:
Adj=(currentSim.Data.AdjMat{IndexTime});%convert 498x498 matrix to EdgeCData ~(1x6065)
Adj=Adj(threshold,threshold);
[j,i,~]=find(tril(Adj));
wd=zeros(1,length(j));
Graph.COMM=Graph.COMM(threshold,threshold);
for k=1:length(j)
    wd(k)=Graph.COMM(i(k),j(k));
end

Adj2=(currentSim.Data.AdjMat{IndexTime});%convert 498x498 matrix to EdgeCData ~(1x6065)
Adj2=Adj2(threshold,threshold);
[j,i,~]=find(tril(Adj2));
wd2=zeros(1,length(j));

for k=1:length(j)
    wd2(k)=Graph.networkThreshold(i(k),j(k));
end

wd3=wd(logical(wd2));
f9=figure;
p9=plot(g);
p9.EdgeCData=wd3;%log10(wd);
p9.MarkerSize=2;
colormap jet
colorbar
labelnode(p9,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
labelnode(p9,highlightElec,{new_electrodes(3).Name{1},new_electrodes(4).Name{1}}); %need to make this better - change 3:4 to a variable

title(['Communicability Analysis Timestamp ' num2str(IndexTime) ' (log10)']);


% CIRCUIT RANK -- measure of recurrent loops (feedback loops)
% based on analyze_network.py
%circuit rank = num edges - num nodes + num connected components
Graph.CircuitRank = numedges(g) - numnodes(g) + sum(conncomp(g));


%Save
network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems

saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Network_Timestamp' num2str(IndexTime)],'jpg');
saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Network_Timestamp' num2str(IndexTime)],'eps');
saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Closeness_Centrality_Timestamp' num2str(IndexTime)],'jpg');
saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Closeness_Centrality_Timestamp' num2str(IndexTime)],'eps');
saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Size_Timestamp' num2str(IndexTime)],'jpg');
saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Size_Timestamp' num2str(IndexTime)],'eps');
saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_&_Closeness_Timestamp' num2str(IndexTime)],'jpg');
saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_&_Closeness_Timestamp' num2str(IndexTime)],'eps');
saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Distribution_Timestamp' num2str(IndexTime)],'jpg');
saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Degree_Distribution_Timestamp' num2str(IndexTime)],'eps');
saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Cluster_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Cluster_Analysis_Timestamp' num2str(IndexTime)],'eps');
saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Participant_Coefficient_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Participant_Coefficient_Analysis_Timestamp' num2str(IndexTime)],'eps');
saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Module_Degree_Z_Score_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Module_Degree_Z_Score_Analysis_Timestamp' num2str(IndexTime)],'eps');
saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Communicability_Analysis_Timestamp' num2str(IndexTime)],'jpg');
saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_Graph_Communicability_Analysis_Timestamp' num2str(IndexTime)],'eps');

end

function save_graph(Graph,network,network_load)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Graph Analysis (Mac)\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_Graph_Analysis_' date],'Graph');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'Graph_Analysis_' num2str(Graph.IndexTime) '_Timestamp_' date],'Graph');
end
end

function save_LDA(LDA_Analysis,network,network_load,sim_loaded)
save_directory='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\LDA Analysis (Mac)\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_LDA_Analysis_' date],'LDA_Analysis');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'LDA_Analysis_' date],'LDA_Analysis','sim_loaded');
end
end
