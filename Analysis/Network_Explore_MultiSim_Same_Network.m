%---------------------------
% 22/04/19
% Network Explore Code v1.1
%
% Author: Alon Loeffler
%
% IMPORTANT:
% Threshold = Degree of greater than 1
% Binarise Threshold = Only low resistence junctions + wires
% Full Graph = all degrees, all resistences
% --------------------------
dbstop if error
%% Load Data

computer=getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        currentPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis';
    case '' %if on linux
        currentPath='/suphys/aloe8475/Documents/CODE/Analysis';
    case 'LAPTOP-S1BV3HR7'
        currentPath='D:\alon_\Research\PhD\CODE\Analysis';
        %case '' %--- Add other computer paths (e.g. Mike)
end
cd(currentPath);
load_data_question=lower(input('Load network data, Analysis Data Only or None? N - None, D - Network Data, A - Analysis Data\n','s'));

if load_data_question=='d'
    clearvars -except load_data_question currentPath
    close all
    %load network data
    [network, network_load, simulations,sim_loaded,numNetworks, explore_network]= load_data(currentPath);
elseif load_data_question=='a'
    clear LDA_Analysis
    close all
    %Load previous LDA analysis data
    networkNum=input(['Which Network # do you want to load? 1 - ' num2str(length(network)) '\n']);
    simNum=input(['Which Simulation # do you want to load? 1 - '  num2str(length(network(networkNum).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    [LDA_Analysis(simNum)] = load_LDA_data(currentPath);
else
    close all
end

%% Choose Simulations
%-------
%Choose Network and Simulation for training or exploring

for currentSimulation=1:length(simulations)
    
if load_data_question~='a'
    if explore_network=='t'
        networkNum=input(['Which Network # do you want to select for Training? 1 - ' num2str(length(network)) '\n']);
        simNum=currentSimulation;%input(['Which Simulation # do you want to select for Training? 1 - '  num2str(network(networkNum).numTrainingSims) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    else
        networkNum=1;%input(['Which Network # do you want to explore? 1 - ' num2str(length(network)) '\n']);
        simNum=currentSimulation;%input(['Which Simulation # do you want to explore? 1 - '  num2str(length(network(networkNum).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    end
end

if size(simulations,1)==1
    currentSim=simulations{simNum};
else
    currentSim=simulations(simNum);
end

% fprintf(['Simulation: ' network(networkNum).Name currentSim.Name ' selected \n\n\n']);


%% Analysis:
i = 1;
while i == 1
    %Choose Analysis to perform
    if explore_network=='t'
        extract_data(network); %extract data for python
        fprintf(['Extracting Data for Python Use...\n\n\n']);
        fprintf(['Data Extracted \n\n\n']);
        analysis_type='e';%lower(input('Which analysis would you like to perform? G - graph, E - Explore Network,L - LDA, N - none \n','s'));
    elseif explore_network=='e' % we don't want to allow LDA if just exploring
        analysis_type='e';%lower(input('Which analysis would you like to perform? G - graph, E - Explore Network, N - none \n','s'));
    end
    if analysis_type=='g'
        %% Graph Analysis
        [Graph,binarise_network]=graph_analysis(network(networkNum),network_load,currentSim,[]);
        %% Save Graph analysis data
        i=i+1;
    elseif analysis_type=='e'
        %% Exploratary analysis of simulation
        
        %Choose a time to Explore Simulation:
        %Find the times where the voltage is on:
        a=currentSim.Data.time(currentSim.Data.VSource1>0.1);
        time1=median(a(1:round(252/3)))-0.005;%252 = total timepoints where voltage >0.5 in all scenarios
        if time1>a(end)
            time1=a(end);
            endTime=1;
        else
            endTime=0;
        end
        time2=time1*5;
        if time2>a(end) & ~endTime
            time2=a(end);
            endTime=1;
        else
            endTime=0;
        end
        time3=time1*10;
        if time3>a(end) & ~endTime
            time3=a(end);
        end
        times=[time1 time2 time3];
        for time=1:length(times)
        [Explore{currentSimulation}(time),threshold{currentSimulation}]=explore_simulation(currentSim,network,network_load,simNum,currentPath,currentSimulation,simulations,times,time);
        end 
        %% Saving Explore
        if currentSimulation==length(simulations)
        save_state=lower(input('Would you like to save the Exploration Analysis? y or n \n','s'));
        if save_state=='y'
            save_explore(Explore,network(networkNum),network_load,currentPath,simNum,threshold,simulations);
            i=i+1; %get out of while loop this loop finishes
        end
        end 
        i=i+1;
    end
    if analysis_type=='g' || analysis_type=='n'
        %% Plotting Graph
        plot_state=lower(input('Would you like to plot Graph Analysis? y or n \n','s'));
        if plot_state=='y'
            Graph=plot_graph(Graph,network(networkNum),network_load,currentSim,sim_loaded,currentPath,binarise_network,simNum);
        end
        i=i+1;
        %% Saving Graph
        save_state=lower(input('Would you like to save the Graph Analysis? y or n \n','s'));
        if save_state=='y'
            save_graph(Graph,network(networkNum),network_load,currentPath);
        end
        i=i+1;
    end
    %% Insert Further Analysis Below
    % -------------------------------
    % -------------------------------
end 
progressBar(currentSimulation,length(simulations));
end 

%--------------------------------------------------------------------------


%% FUNCTIONS

%Loading Functions
function [network, network_load, simulations, sim_loaded, numNetworks, explore_network] = load_data(currentPath)

%% Load Data
%Ask to load Zdenka or Adrian:
network_load='a';%lower(input('Which Network do you want to analyse? Z - Zdenka, A - Adrian \n','s'));

% if strcmp(network_load,'a')
%Get current network - Adrian
[network,sim_loaded, explore_network, numNetworks]=Load_Adrian_Code();
%unpack simulation data into simulation variable
if sim_loaded==1
    if explore_network=='t' %if we have training and testing simulations
        tempSim=network.Simulations{2};
        %         tempSim=num2cell(tempSim);
        %         network.Simulations(2) = [];
        %         network.Simulations=[network.Simulations tempSim];
        
        %number of training + number of testing:
        
        fprintf(['Your Training Simulations are Simulations 1 - ' num2str(network.numTrainingSims) '\n']);
        fprintf(['Your Testing Simulations are Simulations ' num2str(network.numTrainingSims +1) ' - ' num2str(network.numTestingSims) '\n']);
        
        fprintf('\n -------------------------- \nStart Analysis: \n');
        
    end
    for i = 1:length(network.Simulations)
        simulations(i)=network.Simulations(i);
    end
else
    simulations=network.Simulations;
end
cd(currentPath);
% elseif strcmp(network_load,'z')
%Get network - Zdenka:
% D:\alon_\Research\PhD\CODE\Zdenka's Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data
%     network=Load_Zdenka_Code();
%     cd(currentPath);
% end
end

%Exploring Functions:
function [Explore,threshold] = explore_simulation(Sim,network,network_load,simNum,currentPath,currentSimulation,simulations,times,time)

% IMPORTANT:
% Threshold = Degree of greater than 1
% Binarise Threshold = Only low resistence junctions + wires
% Full Graph = all degrees, all resistences

[NodeList.String,NodeList.UserData]=GetNodeList(Sim);
NodeList.Value=1:height(Sim.Electrodes);

%% Timeseries View
%Plot Current
% f=figure;
drainIndex=find(contains(Sim.Electrodes.Name,'Drain'));
if isempty(drainIndex)
    drain_exist=0;
end
sourceIndex=find(contains(Sim.Electrodes.Name,'Source'));
for i = 1:length(sourceIndex)
    source(:,i)=full(Sim.Data.(['ISource' num2str(i)]));
end
for i = 1:length(drainIndex)
    if ismember(['IDrain' num2str(i)], fieldnames(Sim.Data))
        drain(:,i)=full(Sim.Data.(['IDrain' num2str(i)]));
        drain_exist=1;
    end
end
% if drain_exist
%     subplot(1,2,1)
%     plot(source)
%     title('Source')
%     xlabel('Timestamp (0.01sec)')
%     ylabel('Current (A)');
%     subplot(1,2,2)
%     plot(drain)
%     title('Drain');
%     xlabel('Timestamp (0.01sec)')
%     ylabel('Current (A)');
% else
%     plot(source)
%     title('Source')
%     xlabel('Timestamp (0.01sec)')
%     ylabel('Current (A)');
% end

%Plot conductance
% f1=figure;
drainIndex=find(contains(Sim.Electrodes.Name,'Drain'));
if isempty(drainIndex)
    drain_exist=0;
end
sourceIndex=find(contains(Sim.Electrodes.Name,'Source'));
for i = 1:length(sourceIndex)
    source(:,i)=full(Sim.Data.(['ISource' num2str(i)]));
    sourceV(:,i)=full(Sim.Data.(['VSource' num2str(i)]));
end
for i = 1:length(drainIndex)
    if ismember(['IDrain' num2str(i)], fieldnames(Sim.Data))
        drain(:,i)=full(Sim.Data.(['IDrain' num2str(i)]));
        drainV(:,i)=full(Sim.Data.(['VDrain' num2str(i)]));
        drain_exist=1;
    end
end
%     plot(source./sourceV)
%     title('Source')
%     xlabel('Timestamp (0.01sec)')
%     ylabel('Conductance');


IndexTime=single(times(time)*100); %choose the last timestamp 
%input(['What Timestamp do you want to analyse? 1-' num2str(size(Sim.Data,1)) '\n']); %CHOOSE TIMESTAMP

%% Network View
% Function that plots network view of current and resistance
[f2, f3, Adj, NumEl, Explore] = network_view(Sim,IndexTime, NodeList);
% fprintf('Network Analysis Complete \n');

%% Overlay Graph Theory:
% -----------------------------
%Threshold
[Graph, binarise_network]=graph_analysis(network,network_load,Sim,IndexTime);
% fprintf('Graph Analysis Complete \n');
%Threshold graph degree:
if binarise_network=='y'
    threshold_network='t';
else
    threshold_network=lower(input('Do you want to plot the entire Graph or the Thresholded Graph (>1 Degree)? g - entire, t - threshold \n','s'));
end
if threshold_network=='t'
    threshold=Graph.DEG>0; %greater than 0 degree threshold
    Graph.networkThreshold=Graph.AdjMat(threshold,threshold); %applying degree threshold
    G=graph(Graph.networkThreshold);
    node_indices=find(threshold==1); %find nodes with threshold == 1
else
    G=graph(Graph.AdjMat);
end

%% Graph View
% Function that plots graphical view of current, voltage and resistance
if threshold_network=='t'
    [f4, f5, f6, G, Adj, Adj2, Explore, highlightElec, new_electrodes] = graph_view_threshold(Sim,Graph,IndexTime,Explore,G, threshold_network, threshold, drain_exist,node_indices);
else
    [f4, f5, f6, G, Adj, Explore, highlightElec, new_electrodes] = graph_view(Sim,IndexTime,Explore,G, threshold_network,drain_exist,node_indices);
end
%% Graph Theory View
% Function that plots graph theory overlayed on graph view of currents
if threshold_network=='t'
    [f7, f8, f9, f10, f11, f12,f13,f14, Explore,sourceElec, drainElec]= graph_theory_explore_threshold(Sim,G,Adj,Adj2, IndexTime,threshold,threshold_network, Explore, Graph, highlightElec, new_electrodes,node_indices,drain_exist);
%     fprintf('Graph Theory Complete \n');
    
else
    [f7, f8, f9, f10, f11, f12,f13,f14, Explore, sourceElec, drainElec]= graph_theory_explore(Sim,G,Adj,IndexTime,threshold_network, Explore, Graph, highlightElec, new_electrodes,drain_exist);
%     fprintf('Graph Theory Complete \n');
end

 

%% Save
%Save Variables
Explore.IndexTime=IndexTime;
Explore.Name=Sim.Name;
Graph.CircuitRank = numedges(G) - (numnodes(G) - 1);
Explore.GraphTheory=Graph;
Explore.GraphTheory.Definitions={'GE = Global Efficiency','LE = Local Efficiency', 'COMM = Communicability', 'Ci = Community/Cluster Affiliation',...
    'Q = Modularity', 'P = Participation Coefficient', 'MZ = Module Degree z Score', 'AvgPath - Average Path Length', ...
    'GlobalC1ust = number of triangle loops / number of connected triples', 'AvgLocalClust = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)', ...
    'Graph.Clust = a 1xN vector of clustering coefficients per node (where mean(C) = C2)'};
if threshold_network=='t'
    Explore.Thresholded='Yes';
else
    Explore.Thresholded='No';
end
Explore.GraphView.NodeIndices=node_indices;
%% Save Plots
if currentSimulation==length(simulations)
save_explore_plots=lower(input('Would you like to save the plots? y or n \n','s'));
cd(currentPath)
save_directory='..\Data\Figures\Explore Analysis\';
network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems

if save_explore_plots=='y'
    if threshold_network~='t'
        % NOTE: 05/06 - for simulations created after this date, we need to
        % change the file names to Sim.Name, and all this info will be in
        % there already.
%         saveas(f,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Timeseries'],'jpg');
%         print(f,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Timeseries.pdf']);
        saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Conductance_Timeseries' num2str(IndexTime)],'jpg');
        print(f1,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Conductance_Timeseries' num2str(IndexTime) '.pdf']);
        saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
        print(f2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
        print(f3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Currents_Timestamp' num2str(IndexTime)],'jpg');
        print(f4,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
        print(f5,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Voltage_Timestamp' num2str(IndexTime)],'jpg');
        print(f6,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Voltage_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'jpg');
        print(f7,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'jpg');
        print(f8,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Module-zScore_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'jpg');
        print(f9,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Connectivity_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f10,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime)],'jpg');
        print(f10,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f11,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime)],'jpg');
        print(f11,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f12,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime)],'jpg');
        print(f12,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f13,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime)],'jpg');
        print(f13,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f14,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime)],'jpg');
        print(f14,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime) '.pdf']);        
        
    elseif threshold_network=='t'
        saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
        print(f2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
        print(f3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Currents_Timestamp' num2str(IndexTime)],'jpg');
        print(f4,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
        print(f5,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Voltage_Timestamp' num2str(IndexTime)],'jpg');
        print(f6,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Voltage_Timestamp' num2str(IndexTime) '.pdf']);
        saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'jpg');
        print(f7,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'jpg');
        print(f8,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Module-zScore_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'jpg');
        print(f9,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Connectivity_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f10,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime)],'jpg');
        print(f10,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f11,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime)],'jpg');
        print(f11,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f12,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime)],'jpg');
        print(f12,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f13,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime)],'jpg');
        print(f13,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime) '.pdf']);        
        saveas(f14,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime)],'jpg');
        print(f14,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime) '.pdf']);        
           end
end
else
    close all
end

end 

function save_explore(Explore,network,network_load,currentPath,simNum,threshold,Sim)
cd(currentPath);
save_directory='..\Data\Explore Analysis\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_Exploration_Analysis_' date],'Explore');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) '_Sim_' num2str(simNum) '_LastSim_SourceElectrode_' num2str(Explore{simNum}.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore{simNum}.GraphView.ElectrodePosition(2)) '_Exploration_Analysis_ Timestamp_' num2str(Explore{simNum}.IndexTime) '_' date],'Explore','threshold','network','Sim','-v7.3');
end
end

%Graph Functions
function [Graph, binarise_network]=graph_analysis(network,network_load,currentSim,IndexTime)
%% Mac's Analysis: (Graph)
if strcmp(network_load,'z')%Zdenka Code:
    net_mat=network.adj_matrix; %symmetrical matrix
elseif strcmp(network_load,'a') %adrian code
    if isempty(IndexTime)
        IndexTime=input(['What Timestamp do you want to analyse? 1-' num2str(size(currentSim.Data,1)) '\n']); %CHOOSE TIMESTAMP
    end
    %this gives a resistance matrix for the network used for a chosen simulation at a specific timestamp
    binarise_network='y';%input('Do you want to view the network thresholded/binarised ONLY with low Resistance? \n','s');
    if binarise_network=='y' %Binarise so we can use Resistance for graph theory analysis
        
        %what we are doing here is creating a matrix of 0 and 1 (high and
        %low resistence), so that we can plot ONLY those nodes/edges that
        %have low resistence and are relevant.
%         fprintf('Binarising... \n');
        a= full(currentSim.Data.Rmat{IndexTime});
        a(a==5000)=1; %if resistance is low, we make it 1 (on)
        a(a==5000000)=0;  %if it is high we make it 0 (off)
        net_mat=a;
        Graph.binarised='Yes - Using Resistance';
%         fprintf('Binarisation Complete \n');
    else
        fprintf('Binarising... \n');
        net_mat=currentSim.SelLayout.AdjMat; %use standard adjacency matrix
        Graph.binarised='No';
%         fprintf('Binarisation Complete \n');
    end
    %    net_mat=full(simulations(simNum).Data.AdjMat{IndexTime}); %this gives an adj matrix for the network used for a chosen simulation at a specific timestamp
end
logicalAdj=logical(net_mat);
Graph.BinarisedCurrents=currentSim.Data.Currents{IndexTime}(logicalAdj);

%Global efficiency --> 1/characteristic path length, averaged over the whole network. An estimate of how integrated the network is.
% & Distance Matrix
%The distance matrix contains lengths of shortest paths between all pairs of nodes
[Graph.Distance, Graph.GE] = efficiency_bin(net_mat,0);
% fprintf('Global Efficiency Complete \n');

%Topological features will allow us to understand whether it was good or
%bad for classification.

%Local efficiency --> 1/characteristic path length, at each node
[localDistance, Graph.LE]= efficiency_bin(net_mat,1);
% fprintf('Local Efficiency Complete \n');

% Diameter
Graph.Diameter=diameter(net_mat);
% fprintf('Diameter Complete \n');

% Density
%Network Density:
% This is defined, for a given set of nodes, as the number of actual edges
% divided by the number of potential edges.
% I.e., the percentage of possible connections that actually exist.
Graph.Density=density_und(net_mat);
% fprintf('Density Complete \n');

%Betweeness Centrality
Graph.BC=betweenness_bin(net_mat);
% Node betweenness centrality is the fraction of all shortest paths in the network that contain a given node. Nodes with high values of betweenness centrality participate in a large number of shortest paths.
% fprintf('Betweenness Centrality Complete \n');

%Eigenvector Centrality
Graph.EC=eigenvector_centrality_und(full(net_mat));
% Eigenector centrality is a self-referential measure of centrality -- nodes have high eigenvector centrality if they connect to other nodes that have high eigenvector centrality.
% fprintf('Eigenvector Centrality Complete \n');

%Subgraph Centrality
Graph.SC=subgraph_centrality(full(net_mat));
%The subgraph centrality of a node is a weighted sum of closed walks of different lengths in the network starting and ending at the node.
% fprintf('Subgraph Centrality Complete \n');

%communicability --> an estimate of the ease with which each pair of nodes can connect in the network
Graph.COMM = expm(net_mat);
% fprintf('Communicability Complete \n');

%Clustering Coefficient
[Graph.GlobalClust,Graph.AvgLocalClust, Graph.Clust] = clustCoeff(net_mat);
% Graph.Clust2=clustering_coef_bu(net_mat);
% Graph.GlobalC1ust = number of triangle loops / number of connected triples
% Graph.AvgLocalClust = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)
% Graph.Clust = a 1xN vector of clustering coefficients per node (where mean(C) = C2)
% fprintf('Clustering Coefficient Complete \n');


% Graph.Distance=distance_bin(net_mat);
% fprintf('Distance Matrix Complete \n');

%Path Length
Graph.Path = path_length(net_mat);
Graph.CharPath=charpath(Graph.Distance);
%Average Path Length
Graph.AvgPath=mean(Graph.Path);
% fprintf('Path Length Complete \n');

%Connectivity Magnitude Measure of the network:
%Average number of connections divided by maximum number of possible connections 
% In the human brain it is a ratio of ~10^-17 neural connections to possible connections - Kozachkov et al. (2019). How neural circuits achieve and use stable dynamics
Graph.ConnectivityMagnitude=Graph.AvgPath/(length(net_mat)*length(net_mat)); 

%Small World Propensity:
if ~unique(net_mat)==0
Graph.SmallWorldProp=small_world_propensity(net_mat);
end 
% fprintf('Small World Propensity Complete \n');

%modularity --> an estimate of how segregated the network is
[Graph.Ci,Graph.Q] = community_louvain(net_mat,1);
%The Ci term is the module assignment for each node
%The Q term is the 'quality' of the partition --> how modular the network is.
% -- this should tell us how well we can classify
% fprintf('Ci & Q Complete \n');

%Modularity:
if ~unique(net_mat)==0
Graph.Modularity=modularity_und(full(net_mat));
end 
%The optimal community structure is a subdivision of the network into nonoverlapping groups of nodes in a way that maximizes the number of within-group edges, and minimizes the number of between-group edges.
%The modularity is a statistic that quantifies the degree to which the network may be subdivided into such clearly delineated groups
% fprintf('Modularity Complete \n');

%degree --> a count of how many edges are connected to each node
Graph.DEG = degrees_und(net_mat);
% fprintf('Degree Complete \n');

%participation coefficient --> an estimate of how integrative a node is
Graph.P = participation_coef(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'
% fprintf('Participation Coefficient Complete \n');

%module degree z-score --> an estimate of how segregated a node is
Graph.MZ = module_degree_zscore(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'
% fprintf('Module Z Score Complete \n');

%save Adj Matrix
Graph.AdjMat=net_mat;

%save selected time
Graph.IndexTime=IndexTime;
end
function Graph=plot_graph(Graph, network,network_load, currentSim,sim_loaded,currentPath,binarise_network,simNum)
cd(currentPath)
save_directory='..\Data\Figures\Graph Analysis\';

IndexTime=Graph.IndexTime;
%visualise graph network:

%Threshold graph degree:
if binarise_network=='y'
    threshold_choice='t';
else
    lower(input('Do you want to plot the entire Graph or the Thresholded Graph (>1 degree)? g - entire, t - threshold \n','s'));
    threshold_choice=[];
end
if threshold_choice=='t'
    threshold=Graph.DEG>1; %greater than 1 degree threshold - 04/06/19 do I change this to >= 1?
    Graph.networkThreshold=Graph.AdjMat(threshold,threshold); %applying degree threshold
    g=graph(Graph.networkThreshold);
else
    g=graph(Graph.AdjMat);
end
f1=figure;
p1=plot(g);

%find electrodes:
if threshold_choice=='t'
    node_indices=find(threshold==1); %find nodes with threshold == 1
    for i=1:size(currentSim.Electrodes.PosIndex,1)
        if ~isempty(find(node_indices==currentSim.Electrodes.PosIndex(i)))
            new_electrodes(i).PosIndex=find(node_indices==currentSim.Electrodes.PosIndex(i));
            new_electrodes(i).Name=currentSim.Electrodes.Name(i);
        end
    end
else
    for i=1:size(currentSim.Electrodes.PosIndex,1)
        new_electrodes(i).PosIndex=currentSim.Electrodes.PosIndex(i);
        new_electrodes(i).Name=currentSim.Electrodes.Name(i);
    end
end
highlightElec={new_electrodes.PosIndex};
highlightElec=cell2num(highlightElec);
%highlight electrodes on graph
if sim_loaded==1
    highlight(p1,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
    labelnode(p1,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
end

f2=figure;
p2=plot(g);
%% Closeness graph: (Unweighted)
ucc = centrality(g,'closeness');
p2.NodeCData=ucc;
colormap jet
colorbar
title(['Closeness Centrality Scores Timestamp ' num2str(IndexTime)])
if sim_loaded==1
    highlight(p2,highlightElec,'MarkerSize',7); %change simulation number
    if threshold_choice=='t'
        labelnode(p2,[1:size(node_indices,2)],cellstr(num2str(node_indices'))); %label each node with original node number
    end
    labelnode(p2,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
end

%% Node size graph:
f3=figure;
p3=plot(g);
if threshold_choice=='t'
    Graph.DEG_threshold=Graph.DEG(threshold);
    edges = linspace(min(Graph.DEG_threshold),max(Graph.DEG_threshold),7);
    bins = discretize(Graph.DEG_threshold,edges);
else
    edges = linspace(min(Graph.DEG),max(Graph.DEG),7);
    bins = discretize(Graph.DEG,edges);
end
p3.MarkerSize = bins;
p3.NodeColor='r';
title(['Degree Size Timestamp ' num2str(IndexTime)]);
if sim_loaded==1
    highlight(p3,highlightElec,'NodeColor','green'); %change simulation number
    if threshold_choice=='t'
        
        labelnode(p3,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
    end
    labelnode(p3,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
end
text(-5,-6.2,'Min Degrees - 1 (small dot) | Max Degrees - 44 (large dot)');

%% Both combined:
f4=figure;
p4=plot(g);
p4.MarkerSize = bins;
p4.NodeCData=ucc;
colormap jet
colorbar
if threshold_choice=='t'
    labelnode(p4,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
end
labelnode(p4,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
title(['Centrality Score and Degree Size Timestamp ' num2str(IndexTime)]);
text(-5.5,-6.2,'Min Degrees = 1 (small dot) | Max Degrees = 44 (large dot)');

%% Histogram of degree distribution:
f5=figure;
if threshold_choice=='t'
    h1=histogram(Graph.DEG_threshold);
    Graph.avgDEG=mean(Graph.DEG(threshold));
    Graph.stdDEG=std(Graph.DEG(threshold));
else
    h1=histogram(Graph.DEG);
    Graph.avgDEG=mean(Graph.DEG);
    Graph.stdDEG=std(Graph.DEG);
end
title(['Distribution of Connectivity of Nodes Timestamp ' num2str(IndexTime)]);
xlabel('Number of Connections (Degrees)');
ylabel('Frequency');

ylim([0 30]);
text(4,16,['Mean: ' num2str(Graph.avgDEG) ' | SD: ' num2str(Graph.stdDEG)]);

%% Cluster Analysis:
f6=figure;
p6=plot(g);
p6.MarkerSize = 4;
if threshold_choice=='t'
    p6.NodeCData=Graph.Ci(threshold);
    labelnode(p6,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
else
    p6.NodeCData=Graph.Ci;
end
labelnode(p6,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
colormap hsv(13) %change number of colors here if there are more/less than 6 clusters
title(['Cluster Analysis ' num2str(IndexTime)]);

%% Participant Coefficient Analysis:
f7=figure;
p7=plot(g);
if threshold_choice=='t'
    p7.NodeCData=Graph.Ci(threshold);
    p_ranks=Graph.P(threshold);
else
    p7.NodeCData=Graph.Ci;
    p_ranks=Graph.P;
end
edges2 = linspace(min(p_ranks),max(p_ranks),7);
bins2 = discretize(p_ranks,edges2);
p7.MarkerSize=bins2;
colormap hsv(6)
if threshold_choice=='t'
    labelnode(p7,[1:size(node_indices,2)],cellstr(num2str(node_indices'))); %label each node with original node number
end
labelnode(p7,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
if threshold_choice=='t'
    text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P(threshold))) ' (large dot)']);
else
    text(-5.5,-6.2,['Min P Coeff = 0 (small dot) | Max P Coeff = ' num2str(max(Graph.P)) ' (large dot)']);
end
title(['Participant Coefficient Analysis Timestamp ' num2str(IndexTime)]);

%% Module Degree Z-Score:
f8=figure;
p8=plot(g);
if threshold_choice=='t'
    
    p8.NodeCData=Graph.Ci(threshold);
    mod_ranks=Graph.MZ(threshold);
else
    
    p8.NodeCData=Graph.Ci;
    mod_ranks=Graph.MZ;
end
edges3 = linspace(min(mod_ranks),max(mod_ranks),7);
bins3 = discretize(mod_ranks,edges3);
p8.MarkerSize=bins3;
colormap hsv(6)
if threshold_choice=='t'
    labelnode(p8,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
end
labelnode(p8,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable
if threshold_choice=='t'
    text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ(threshold))) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ(threshold))) ' (large dot)']);
else
    text(-6,-6.2,['Min MZ Coeff = ' num2str(min(Graph.MZ)) ' (small dot) | Max MZ Coeff = ' num2str(max(Graph.MZ)) ' (large dot)']);
end
title(['Module Degree z-Score Analysis Timestamp ' num2str(IndexTime)]);


%% Communicability at different times:
Adj=(currentSim.Data.AdjMat{IndexTime});%convert 498x498 matrix to EdgeCData ~(1x6065)
if threshold_choice=='t'
    Adj=Adj(threshold,threshold);
    Graph.COMM=Graph.COMM(threshold,threshold);
end
% extract lower triangular part of Adjacency matrix of Graph.COMM
[j,i,~]=find(tril(Adj));
wd=zeros(1,length(j));

for k=1:length(j)
    wd(k)=Graph.COMM(i(k),j(k));
end

% extract lower triangular part of Adjacency matrix of network
Adj2=(currentSim.Data.AdjMat{IndexTime});%convert 498x498 matrix to EdgeCData ~(1x6065)
if threshold_choice=='t'
    
    Adj2=Adj2(threshold,threshold);
end
[j,i,~]=find(tril(Adj2));
wd2=zeros(1,length(j));

for k=1:length(j)
    if threshold_choice=='t'
        wd2(k)=Graph.networkThreshold(i(k),j(k));
    else
        wd2(k)=Graph.AdjMat(i(k),j(k));
    end
end

%Find Graph.COMM in network
wd3=wd(logical(wd2));

f9=figure;
p9=plot(g);
p9.EdgeCData=wd3;%log10(wd);
p9.MarkerSize=2;
colormap jet
colorbar
if threshold_choice=='t'
    labelnode(p9,[1:size(node_indices,2)],cellstr(num2str(node_indices')));  %label each node with original node number
end
labelnode(p9,highlightElec,[new_electrodes(:).Name]); %need to make this better - change 3:4 to a variable

title(['Communicability Analysis Timestamp ' num2str(IndexTime) ' (log10)']);


%% CIRCUIT RANK -- measure of recurrent loops (feedback loops)
% based on analyze_network.py
%circuit rank = num edges - num nodes + num connected components
Graph.CircuitRank = numedges(g) - (numnodes(g) - 1);
if threshold_choice=='t'
    Graph.Indices=node_indices;
end

%% Save
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
function save_graph(Graph,network,network_load,currentPath)
cd(currentPath);
save_directory='..\Data\Graph Analysis (Mac)\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_Graph_Analysis_' date],'Graph');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'Graph_Analysis_' num2str(Graph.IndexTime) '_Timestamp_' date],'Graph');
end
end

