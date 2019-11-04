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

set(0,'DefaultFigureVisible','off')

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
load_data_question=lower(input('Load network data, Analysis Data Only or None? N - None, D - Network Data, A - Analysis Data, L - All Network Data for Loop \n','s'));

if load_data_question=='d' || load_data_question=='l'
    if load_data_question=='l'
        loop=1;
    else
        loop=0;
    end
    clearvars -except load_data_question currentPath loop
    close all
    %load network data
    [network, network_load, simulations,sim_loaded,numNetworks, explore_network]= load_data(currentPath,loop);
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
if load_data_question~='a'
    if explore_network=='t'
        networkNum=input(['Which Network # do you want to select for Training? 1 - ' num2str(length(network)) '\n']);
        simNum=input(['Which Simulation # do you want to select for Training? 1 - '  num2str(network(networkNum).numTrainingSims) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    else
        networkState=input(['Which Network # do you want to explore?  1  -  ' num2str(length(network)) '  OR  0  to loop through all \n']);
        networkNum=networkState;
        if networkState==0
            simNum=1;
        else
            simNum=input(['Which Simulation # do you want to explore? | 1 | - | '  num2str(length(network(networkNum).Simulations)) ' | \n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
        end
    end
end

%% Loop through networks:
if networkNum==0 %if we want to loop through all networks
    NetMulti=[];
    loop=1;
    for networkCount=1:length(network) %loop though the networks:
        currentSim=network(networkCount).Simulations{simNum};
        % ANALYSIS: -----------------------------------------------------
        i = 1;
        while i == 1
            %% Exploratary analysis of simulation
            Explore{networkCount}=explore_simulation_loop(currentSim,network(networkCount),network_load,simNum,currentPath,loop,networkCount);
            i=i+1;
            %% Saving Explore
            if networkCount==length(network)
                save_state=lower(input('Would you like to save the Exploration Analysis? y or n \n','s'));
                if save_state=='y'
                    save_explore(Explore,network(networkCount),network_load,currentPath,simNum,loop,currentSim);
                    i=i+1; %get out of while loop this loop finishes
                end
            end
            %% Insert Further Analysis Below
            % -------------------------------
            % -------------------------------
        end
    end
else %if we're not looping
    if sim_loaded
%         if size(simulations,1)==1
%             currentSim=simulations{simNum};
%         else
            currentSim=simulations(simNum);
%         end
    else
        currentSim=network(networkNum).Simulations{simNum};
    end
    fprintf(['Simulation: ' network(networkNum).Name  ' selected \n\n\n']); %currentSim.Name
    
    
    %% No Looping through networks
    % Analysis -----------------------------------------------
    i = 1;
    while i == 1
        %Choose Analysis to perform
        if explore_network=='t'
            extract_data(network); %extract data for python
            fprintf(['Extracting Data for Python Use...\n\n\n']);
            fprintf(['Data Extracted \n\n\n']);
            analysis_type=lower(input('Which analysis would you like to perform? G - graph, E - Explore Network,L - LDA, N - none \n','s'));
        elseif explore_network=='e' % we don't want to allow LDA if just exploring
            analysis_type=lower(input('Which analysis would you like to perform? G - graph, E - Explore Network, N - none \n','s'));
        end
        if analysis_type=='g'
            %% Graph Analysis
            [Graph,binarise_network]=graph_analysis(network(networkNum),network_load,currentSim,[],loop);
            %% Save Graph analysis data
            i=i+1;
        elseif analysis_type=='e'
            %% Exploratary analysis of simulation
            Explore=explore_simulation(currentSim,network,network_load,simNum,currentPath,loop);
            %% Saving Explore
            save_state=lower(input('Would you like to save the Exploration Analysis? y or n \n','s'));
            if save_state=='y'
                save_explore(Explore,network(networkNum),network_load,currentPath,simNum,loop,currentSim);
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
                plot_LDA(LDA_Analysis(simNum),simNum,network(networkNum).Name,currentPath,simulations);
                i=i+1;
            end
            
            %% Apply LDA Training to testing data (different simulation)
            if analysis_type~='n'
                apply_LDA=lower(input('Would you like to apply the loaded LDA analysis to another simulation? y or n \n','s'));
                if apply_LDA=='y'
                    [LDA_Analysis, simulationChoice]=lda_apply_func(numNetworks,network,LDA_Analysis,simNum,simulations,currentPath);
                end
                i=i+1; %get out of while loop when this loop finishes
            elseif analysis_type~='l' &&  analysis_type~='g' && analysis_type~='n' && analysis_type~='e'
                fprintf('Please type either G, L or N only \n');
            end
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
end

%--------------------------------------------------------------------------


%% FUNCTIONS

%Loading Functions
function [network, network_load, simulations, sim_loaded, numNetworks, explore_network] = load_data(currentPath,loop)

%% Load Data
%Ask to load Zdenka or Adrian:
network_load=lower(input('Which Network do you want to analyse? Z - Zdenka, A - Adrian \n','s'));

if ~loop
    if strcmp(network_load,'a')
        
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
        % D:\alon_\Research\PhD\CODE\Zdenka Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data
        %     network=Load_Zdenka_Code();
        %     cd(currentPath);
        % end
        
    elseif strcmp(network_load,'z')
        % Get network - Zdenka:
        % D:\alon_\Research\PhD\CODE\Zdenka Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data
        [network,sim_loaded, explore_network, numNetworks]=Load_Zdenka_Code();
        simulations=network.Simulations;
        cd(currentPath);
        % end
    end
    
else
    % Load all data from specified folder:
    currMultiNet=input('Which size Network would you like to load iterations from? 100, 500, 1000 or 2000? \n');
    computer=getenv('computername');
    switch computer
        case 'W4PT80T2'
            dataPath=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Adrian Networks\' num2str(currMultiNet) 'nw Alternate NWs\'];
        case ''
            dataPath=['/suphys/aloe8475/Documents/CODE/Adrian''s Code/NETWORK_sims_2/Saved Networks/' num2str(currMultiNet) 'nw Alternate NWs/'];
        case 'LAPTOP-S1BV3HR7'
            dataPath=['D:\alon_\Research\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\' num2str(currMultiNet) 'nw Alternate NWs\'];
    end
    fprintf(['Loading Multi Network Data from Network size ' num2str(currMultiNet) '\n']);
    cd(dataPath); %go to chosen data folder
    dinfo = dir('Net*.mat'); %choose files starting with 'Net' in folder.
    for K = 1 : length(dinfo) %Load all files in selected folder
        thisfile = dinfo(K).name;
        destfile = fullfile(dataPath, thisfile);
        temp = load(thisfile); %load file into network variable
        if strcmp(network_load,'a') %if adrian data structure
            network(K)=temp.SelNet;
        elseif strcmp(network_load,'z') %if zdenka data structure
            network(K)=temp.SelNet;
        end
        simulations(K)=network(K).Simulations;
        clear temp
    end
    numNetworks=length(network); %number of networks
    sim_loaded=1; %yes we are loading simulations
    explore_network='e'; %explore
end


end
function [LDA_Analysis] = load_LDA_data(currentPath)
cd(currentPath);
cd('..\Data\LDA Analysis (Mac)');
waitfor(msgbox('Select the LDA Analysis saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the LDA Analysis saved data');
f=fullfile(PathName,FileName);
load(f);
cd('..\..\CODE\Analysis');
end

%Exploring Functions:
function Explore=explore_simulation_loop(Sim,network,network_load,simNum,currentPath,loop,networkNum)
% IMPORTANT:
% Threshold = Degree of greater than 1
% Binarise Threshold = Only low resistence junctions + wires
% Full Graph = all degrees, all resistences
[NodeList.String,NodeList.UserData]=GetNodeList(Sim);
NodeList.Value=1:height(Sim.Electrodes);
fprintf('Starting Network Analysis... \n');
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
% plot(source./sourceV)
% title('Source')
% xlabel('Timestamp (0.01sec)')
% ylabel('Conductance');

%Choose a time to Explore Simulation:
%Always choose the last available time stamp:
IndexTime=size(Sim.Data,1);%input(['What Timestamp do you want to analyse? 1-' num2str(size(Sim.Data,1)) '\n']); %CHOOSE TIMESTAMP

%% Network View
% Function that plots network view of current and resistance
% [f2, f3, Adj, NumEl, Explore] = network_view(Sim,IndexTime, NodeList,network_load);
[Adj, NumEl, Explore] = Multinetwork_view(Sim,IndexTime, NodeList,network_load);

fprintf('Network View Analysis Complete \n');
close all
%% Overlay Graph Theory:
% -----------------------------
%Threshold
[Graph, binarise_network]=graph_analysis(network,network_load,Sim,IndexTime,loop);
fprintf('Graph Analysis Complete \n');
%Threshold graph degree:
if binarise_network=='y'
    threshold_network='t';
else
    threshold_network='g';%lower(input('Do you want to plot the entire Graph or the Thresholded Graph (>1 Degree)? g - entire, t - threshold \n','s'));
end
if threshold_network=='t'
    threshold=Graph.DEG>0; %greater than 1 degree threshold
    Graph.networkThreshold=Graph.AdjMat(threshold,threshold); %applying degree threshold
    G=graph(Graph.networkThreshold);
    node_indices=find(threshold==1); %find nodes with threshold == 1
else
    G=graph(Graph.AdjMat);
end

%% Graph View
% Function that plots graphical view of current, voltage and resistance
if threshold_network=='t'
    [f4, f5, f6, G, Adj, Adj2, Explore, highlightElec, new_electrodes] = graph_view_threshold(Sim,Graph,IndexTime,Explore,G, threshold_network, threshold, drain_exist);
else
    [f4, f5, f6, G, Adj, Explore, highlightElec, new_electrodes] = Multigraph_view(Sim,IndexTime,Explore,G, threshold_network,drain_exist,network_load);
end
close all
%% Graph Theory View
% Function that plots graph theory overlayed on graph view of currents
if threshold_network=='t'
    [f7, f8, f9, f10, f11, f12,f13, Explore,sourceElec, drainElec]= graph_theory_explore_threshold(Sim,G,Adj,Adj2, IndexTime,threshold,threshold_network, Explore, Graph, highlightElec, new_electrodes,node_indices,drain_exist,network_load);
    fprintf('Graph Theory Complete \n');
    
else
    [f7, f8, f9, f10, f11, f12,f13, Explore, sourceElec, drainElec]= Multigraph_theory_explore(Sim,G,Adj,IndexTime,threshold_network, Explore, Graph, highlightElec, new_electrodes,drain_exist,[],network_load);
    fprintf('Graph Theory Complete \n');
end
close all

fG=figure;
pG=plot(Explore.GraphView.G);
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

fprintf(['\n \n Network ' num2str(networkNum) ' Finished \n \n']);
end
function Explore = explore_simulation(Sim,network,network_load,simNum,currentPath,loop)

% IMPORTANT:
% Threshold = Degree of greater than 0
% Binarise Threshold = Only low resistence junctions + wires
% Full Graph = all degrees, all resistences

[NodeList.String,NodeList.UserData]=GetNodeList(Sim,network_load);
NodeList.Value=1:length([Sim.Electrodes.PosIndex]);

%% Timeseries View
%Plot Current
f=figure;
drainIndex=find(contains(NodeList.String,'Drain'));
if isempty(drainIndex)
    drain_exist=0;
end
sourceIndex=find(contains(NodeList.String,'Source'));
for i = 1:length(sourceIndex)
    if ismember(['ISource' num2str(i)],fieldnames(Sim.Data))
    source(:,i)=full(Sim.Data.(['ISource' num2str(i)]));
    sourceV(:,i)=full(Sim.Data.(['VSource' num2str(i)]));
    source_exist=1;
    end 
end
for i = 1:length(drainIndex)
    if ismember(['IDrain' num2str(i)], fieldnames(Sim.Data))
        drain(:,i)=full(Sim.Data.(['IDrain' num2str(i)]));
        drainV(:,i)=full(Sim.Data.(['VDrain' num2str(i)]));
        drain_exist=1;
    end
end
if drain_exist
    subplot(1,2,1)
    plot(-source)
    title('Source')
    xlabel('Timestamp (0.01sec)')
    ylabel('Current (A)');
    subplot(1,2,2)
    plot(-drain)
    title('Drain');
    xlabel('Timestamp (0.01sec)')
    ylabel('Current (A)');
else
    plot(-source)
    title('Source')
    xlabel('Timestamp (0.01sec)')
    ylabel('Current (A)');
end

%Plot conductance
f1=figure;
plot(-source./sourceV)
title('Source')
xlabel('Timestamp (0.01sec)')
ylabel('Conductance');

%Choose a time to Explore Simulation:
IndexTime=input(['What Timestamp do you want to analyse? 1-' num2str(length(Sim.Time)) '\n']); %CHOOSE TIMESTAMP

%% Network View
% Function that plots network view of current and resistance
% [f2, f3, Adj, NumEl, Explore] = network_view(Sim,IndexTime, NodeList,network_load);
% NEED TO COMMENT THIS OUT UNTIL CAN FIX CY AND CX

Explore=[];

fprintf('Network Analysis Complete \n');
%% Overlay Graph Theory:
% -----------------------------
%Threshold
[Graph, binarise_network]=graph_analysis(network,network_load,Sim,IndexTime,loop);
fprintf('Graph Analysis Complete \n');
%Threshold graph degree:
if binarise_network=='y'
    threshold_network='t';
else
    threshold_network=lower(input('Do you want to plot the entire Graph or the Thresholded Graph (>1 Degree)? g - entire, t - threshold \n','s'));
end
if threshold_network=='t'
    threshold=Graph.DEG>0; %greater than 1 degree threshold
    Graph.networkThreshold=Graph.AdjMat(threshold,threshold); %applying degree threshold
    G=graph(Graph.networkThreshold);
    node_indices=find(threshold==1); %find nodes with threshold == 1
else
    G=graph(Graph.AdjMat);
end

%% Graph View
% Function that plots graphical view of current, voltage and resistance
if threshold_network=='t'
    [f4, f5, f6, G, Adj, Adj2, Explore, highlightElec, new_electrodes] = graph_view_threshold(Sim,Graph,IndexTime,Explore,G, threshold_network, threshold, drain_exist,network_load,node_indices);
else
    [f4, f5, f6, G, Adj, Explore, highlightElec, new_electrodes] = graph_view(Sim,IndexTime,Explore,G, threshold_network,drain_exist,network_load);
end

%% Graph Theory View
% Function that plots graph theory overlayed on graph view of currents
if threshold_network=='t'
    [f7, f8, f9, f10, f11, f12,f13, f14, Explore,sourceElec, drainElec]= graph_theory_explore_threshold(Sim,G,Adj,Adj2, IndexTime,threshold,threshold_network, Explore, Graph, highlightElec, new_electrodes,node_indices,drain_exist,source_exist,network_load);
    fprintf('Graph Theory Complete \n');
    
else
    [f7, f8, f9, f10, f11, f12,f13, f14, Explore, sourceElec, drainElec]= graph_theory_explore(Sim,G,Adj,IndexTime,threshold_network, Explore, Graph, highlightElec, new_electrodes,drain_exist,source_exist,network_load);
    fprintf('Graph Theory Complete \n');
end


fG=figure;
pG=plot(G);
pG.NodeLabel=[];


%Biograph view
% h = view(biograph(Adj,[],'ShowArrows','off'));
% set(h.Nodes(path),'Color',[1 0.4 0.4])
% fowEdges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
% revEdges = getedgesbynodeid(h,get(h.Nodes(fliplr(path)),'ID'));
% edges = [fowEdges;revEdges];
% set(edges,'LineColor',[1 0 0])
% set(edges,'LineWidth',1.5)

%% Searching Algorithms

if threshold_network=='t' %only conduct search if we thresholded the network - otherwise too complex
    %     T1 = bfsearch(G,sourceElec,'allevents'); %Breadth-First Search
    %     T2 = dfsearch(G, sourceElec, 'allevents', 'Restart', true); %Depth-First Search
    %figure;
    %visualize_search(G,T1) %Visual search step by step
    % visualize_search(G,T2) %Visual search step by step
    
    %plot searching algorithms
    %bfsearch
    %     fs1=figure;p = plot(G,'Layout','layered');
    %     events = {'edgetonew','edgetofinished','startnode'};
    %     T = bfsearch(G,sourceElec,events,'Restart',true);
    %     highlight(p, 'Edges', T.EdgeIndex(T.Event == 'edgetofinished'), 'EdgeColor', 'k')
    %     highlight(p, 'Edges', T.EdgeIndex(T.Event == 'edgetonew'), 'EdgeColor', 'r')
    %     highlight(p,T.Node(~isnan(T.Node)),'NodeColor','g')
    %     if drain_exist
    %     %Overlay shortest path:
    %     [dist,path,pred]=graphshortestpath(Adj2,sourceElec,drainElec,'Directed','false');
    %     highlight(p,path,'EdgeColor','cyan','LineWidth',6,'LineStyle','-');
    %     title('Layered Graph Breadth-First Search overlayed w Shortest Path');
    %         fprintf('Shortest Path Complete \n');
    %     end
    %Outputs:
    
    %'discovernode' (default)-A new node has been discovered.
    %'finishnode'- All outgoing edges from the node have been visited.
    %'startnode'- This flag indicates the starting node in the search. If 'Restart' is true, then 'startnode' flags the starting node each time the search restarts.
    %'edgetonew'-Edge connects to an undiscovered node.
    %'edgetodiscovered'	-Edge connects to a previously discovered node.
    %'edgetofinished'- Edge connects to a finished node.
    
end

%% Save
%Save Variables
Explore.IndexTime=IndexTime;
if network_load=='a'
    Explore.Name=Sim.Name;
end

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
%% Save Plots

save_explore_plots=lower(input('Would you like to save the plots? y or n \n','s'));
cd(currentPath)
save_directory='..\Data\Figures\Explore Analysis\';
if length(network)== 1
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
else
    for q = 1:length(network)
        network(q).Name(regexp(network(q).Name,'[/:]'))=[];
    end
end
if save_explore_plots=='y'
    if threshold_network~='t'
        % NOTE: 05/06 - for simulations created after this date, we need to
        % change the file names to Sim.Name, and all this info will be in
        % there already.
        %         saveas(f,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Timeseries'],'jpg');
        %         print(f,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Timeseries.pdf']);
        print(fG,'-painters','-dpng','-r300',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_ExploreGraph.png']);
        print(fG,'-painters','-dpdf','-bestfit','-r300',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_ExploreGraph.pdf']);
        
        saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Conductance_Timeseries' num2str(IndexTime)],'jpg');
        print(f1,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Conductance_Timeseries' num2str(IndexTime) '.pdf']);
%         saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
%         print(f2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
% %         saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
%         print(f3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
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
        %         saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
%         print(f2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
%         saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
%         print(f3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
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
end
function save_explore(Explore,network,network_load,currentPath,simNum,loop,currentSim)
cd(currentPath);
% save_directory='..\Data\Explore Analysis\';
if loop==0
    save_directory=['..\Data\Explore Analysis\Continuous DC\'];
elseif loop ==1
    if network_load=='z'
        save_directory=['..\Data\Explore Analysis\MultiNetwork Analysis\' num2str(network.NetworkSettings.Number) 'nw Alternate NWs\'];
    else
        save_directory=['..\Data\Explore Analysis\MultiNetwork Analysis\' num2str(network.NetworkSettings.Number) 'nw Alternate NWs\'];
    end 
end
if strcmp(network_load,'z')%Zdenka Code:
        network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Zdenka_' num2str(network.Name) 'Length_' num2str(network.NetworkSettings.Length) '_Disp_' num2str(network.NetworkSettings.Disp) '_Sim_' num2str(simNum) '_Source_' num2str(Explore{1}.GraphView.ElectrodePosition(1)) '_Drain_' num2str(Explore{1}.GraphView.ElectrodePosition(2)) '_Explore_Timestamp_' num2str(Explore{1}.IndexTime)],'Explore');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    if loop
        save([save_directory 'Adrian_NW_Sx' num2str(network.NetworkSettings.Number) '_NumNetworkIterations_' num2str(length(Explore)) '_' date],'Explore');
    else
        save([save_directory 'Adrian_' num2str(network.Name) '_Length_' num2str(network.NetworkSettings.Length) '_Disp_' num2str(network.NetworkSettings.Disp) '_Sim_' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Exploration_Analysis_Timestamp_' num2str(Explore.IndexTime)],'Explore');
    end
end
end

%LDA Functions
function LDA_Analysis=lda_analysis(currentSim,network,network_load,simNum)
%% LDA Analysis
%THIS IS WHERE YOU CHANGE THE VARIABLES FOR LDA
drain1=[full(currentSim.Data.IDrain1)]; %full(simulations(2).Data.IDrain1)];
if exist('currentSim.Data.IDrain2','var')
    drain2=[full(currentSim.Data.IDrain2)]; %full(simulations(2).Data.IDrain2)];
    exists =1;
else
    exists =0;
end
%Current source
source1=[full(currentSim.Data.ISource1)]; %full(simulations(2).Data.ISource1)];
source2=[full(currentSim.Data.ISource2)]; %full(simulations(2).Data.ISource2)];
%Voltage source
source1Voltage=full(currentSim.Data.VSource1);
source2Voltage=full(currentSim.Data.VSource2);
Input = [source1 source2];
if exists
    Output = [drain1 drain2]; %OUTPUT is column (variable) x row (observations)
else
    Output=drain1;
end
Target = source1Voltage(:,1)>0.001;% uncomment for current: [source1(:,1)> 1.0e-04 *0.001 & source2<0];%TARGET is the classifier we expect

%visualise drain1 and drain2
figure
plot(Output);

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
if exists
    LDA_Analysis(simNum).drain2=drain2;
end
LDA_Analysis(simNum).source1=source1;
LDA_Analysis(simNum).source2=source2;
LDA_Analysis(simNum).Output=Output;
LDA_Analysis(simNum).Input=Input;
LDA_Analysis(simNum).Target=Target;

LDA_Analysis(simNum).normalisedOutput=LDA_normalise_drain(Output);
LDA_Analysis(simNum).normalisedInput=LDA_normalise_source(Input);
LDA_Analysis(simNum).normalisedW=LDA(LDA_Analysis(simNum).normalisedOutput,Target);
LDA_Analysis(simNum).normalisedL=[ones(size(LDA_Analysis(simNum).normalisedOutput,1),1) LDA_Analysis(simNum).normalisedOutput] * LDA_Analysis(simNum).normalisedW';
LDA_Analysis(simNum).normalisedP=exp(LDA_Analysis(simNum).normalisedL) ./ repmat(sum(exp(LDA_Analysis(simNum).normalisedL),2),[1 2]);

LDA_Analysis(simNum).TypeOfData='Training';


% %% Support Vector Machine Analysis:
% rng default
%
% SVMModel = fitcsvm(Output,Target,'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus'));
% sv = SVMModel.SupportVectors;
% figure
% gscatter(Output(:,1),Output(:,2),Target)
% hold on
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% legend('Source 1','Source 2','Support Vector','Location','West')
% hold off

%Cross Validate:
% CVSVMModel = crossval(SVMModel);
% classLoss = kfoldLoss(CVSVMModel); %Generalization Rate


if exists
    clear drain1 drain2 source1 source2 Input Target Output
else
    clear drain1 source1 source2 Input Target Output
end
%% Saving LDA
save_state=lower(input('Would you like to save the LDA Analysis? y or n \n','s'));
if save_state=='y'
    save_LDA(LDA_Analysis(simNum),network(networkNum),network_load,sim_loaded,currentPath);
end
end
function [LDA_Analysis, simulationChoice]=lda_apply_func(numNetworks,network,LDA_Analysis,simNum,simulations,currentPath)
if numNetworks>1 %if we have two networks, offer to test second network
    networkNum2=input(['Which Network # do you want to select your Simulation from ? 1 - ' num2str(length(network)) '\n']);
else
    networkNum2=1;
end
if numNetworks>1 %if we have more than 1 simulation, they can input 1 as an option
    simulationChoice=input(['Which Simulation # do you want to apply LDA to? 1 - '  num2str(length(network(networkNum2).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
else % otherwise, can only input 2 or higher
    while 1
        simulationChoice=input(['Which Simulation # do you want to apply LDA to? ' num2str(network(networkNum2).numTrainingSims+1) '- '  num2str(network(networkNum2).numTestingSims) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
        if simulationChoice >1 && simulationChoice <=length(network(networkNum2).Simulations)
            break;
        elseif simulationChoice ==1
            fprintf('Cannot Train and Test the same Simulation, please choose again \n');
        else
            fprintf(['Please Choose a number between '  num2str(network(networkNum2).numTrainingSims+1) ' and ' num2str(network(networkNum2).numTestingSims)]);
        end
    end
    if exist('simulations{simulationChoice}.Data.IDrain2','var')
        LDA_Analysis(simulationChoice).Output=[full(simulations{simulationChoice}.Data.IDrain1) full(simulations{simulationChoice}.Data.IDrain2)];
    else
        LDA_Analysis(simulationChoice).Output=[full(simulations{simulationChoice}.Data.IDrain1)];
    end
    LDA_Analysis(simulationChoice).Input=[full(simulations{simulationChoice}.Data.ISource1) full(simulations{simulationChoice}.Data.ISource2)];
    LDA_Analysis(simulationChoice).Target=[full(simulations{simulationChoice}.Data.VSource1 >0.001)];
    LDA_Analysis(simulationChoice).TypeOfData='Testing';
    [LDA_Analysis(simulationChoice).appliedP, LDA_Analysis(simulationChoice).appliedL, LDA_Analysis(simulationChoice).normalisedOutput, LDA_Analysis(simulationChoice).normalisedInput]=LDA_Apply(LDA_Analysis(simNum).normalisedW,LDA_Analysis(simulationChoice).Output, LDA_Analysis(simulationChoice).Input);
    plot_state2=lower(input('Would you like to plot the applied LDA results? y or n \n','s'));
    if plot_state2=='y'
        SelSims=simulations{simNum};
        plot_LDA_Applied(LDA_Analysis(simulationChoice),simNum,simulationChoice,network(networkNum2).Name,currentPath,SelSims);
    end
end
end
function plot_LDA(LDA_Analysis, simNum, networkName,currentPath,simulations)
cd(currentPath);
SelSims=simulations{simNum};
save_directory='..\Data\Figures\LDA\LDA Training\';
%Plot Electrode Locations:
fElectrodes=figure;
Sim=simulations{simNum};
AdjMat=Sim.SelLayout.AdjMat;
G=graph(AdjMat);
p1=plot(G);
p1.NodeLabel={};
highlightElec=Sim.Electrodes.PosIndex;
highlight(p1,highlightElec,'NodeColor','green','MarkerSize',5); %change simulation number
labelnode(p1,highlightElec,[Sim.Electrodes.Name]); %need to make this better - change 3:4 to a variable

title('Training Electrodes (Graph View)');

%plot inputs
fin=figure;
plot(LDA_Analysis.Input)
title('Source Electrodes');
ylabel('Current (A)');
xlabel('Timestamp (0.01sec)');
%plot outputs
fout=figure;
plot(LDA_Analysis.Output)
title('Drain Electrodes');
ylabel('Current (A)');
xlabel('Timestamp (0.01sec)');
%plot LDA
LDAf=figure;
sgtitle('LDA Classification Training');
subplot(4,1,1)
plot(LDA_Analysis.Input)
title('Source Electrodes');
ylabel('Current (A)');
xlabel('Timestamp (0.01sec)');
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

%Plot LDA Scatter + Seperator:



%Save
%Note; change the date of the simulation in the name if using a different
%set of simulations

save_status=lower(input('Would you like to save the LDA Plots? y or n \n','s'));
if save_status=='y'
    SelSims.Name=strrep(SelSims.Name,'.','');
    saveas(fElectrodes,[save_directory num2str(networkName) SelSims.Name  '_Sec_Electrode_Position'],'jpg');
    saveas(fElectrodes,[save_directory num2str(networkName) SelSims.Name  '_Sec_Electrode_Position'],'eps');
    saveas(fout,[save_directory num2str(networkName) SelSims.Name  '_Sec_Training_Inputs'],'jpg');
    saveas(fout,[save_directory num2str(networkName) SelSims.Name  '_Sec_Training_Inputs'],'eps');
    saveas(fin,[save_directory num2str(networkName) SelSims.Name   '_Sec_Training_Outputs'],'jpg');
    saveas(fin,[save_directory num2str(networkName) SelSims.Name  '_Sec_Training_Outputs'],'eps');
    saveas(LDAf,[save_directory num2str(networkName) SelSims.Name  '_Sec_LDA_Classification_Training'],'jpg');
    saveas(LDAf,[save_directory num2str(networkName) SelSims.Name   '_Sec_LDA_Classification_Training'],'eps');
    saveas(LDAff,[save_directory num2str(networkName) SelSims.Name   '_Sec_logLDA_Classification_Training'],'jpg');
    saveas(LDAff,[save_directory num2str(networkName) SelSims.Name  '_Sec_logLDA_Classification_Training'],'eps');
    saveas(LDAnormf,[save_directory num2str(networkName) SelSims.Name  '_Sec_norm_LDA_Classification_Training'],'jpg');
    saveas(LDAnormf,[save_directory num2str(networkName) SelSims.Name  '_Sec_norm_LDA_Classification_Training'],'eps');
end
end
function plot_LDA_Applied(LDA_Analysis,simNum,simulationChoice,networkName,currentPath,SelSims)
cd(currentPath)
save_directory='..\Data\Figures\LDA\LDA Testing\';
%plot inputs
fin=figure;
plot(LDA_Analysis.Input)
title('Source Electrodes');
ylabel('Current (A)');
xlabel('Timestamp (0.01sec)');
%plot outputs
fout=figure;
plot(LDA_Analysis.Output)
title('Drain Electrodes');
ylabel('Current (A)');
xlabel('Timestamp (0.01sec)');

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
save_status=lower(input('Would you like to save the Applied LDA Plots? y or n \n','s'));


if save_status=='y'
    SelSims.Name=strrep(SelSims.Name,'.','');
    
    saveas(fout,[save_directory num2str(networkName) SelSims.Name '_Sec_Testing_Inputs'],'jpg'); % CHANGE DATE OF SIMULATION
    saveas(fout,[save_directory num2str(networkName) SelSims.Name '_Sec_Testing_Inputs'],'eps');
    saveas(fin,[save_directory num2str(networkName) SelSims.Name  '_Sec_Testing_Outputs'],'jpg'); % CHANGE DATE OF SIMULATION
    saveas(fin,[save_directory num2str(networkName) SelSims.Name  '_Sec_Testing_Outputs'],'eps');
    
    saveas(appliedF,[save_directory num2str(networkName) SelSims.Name '_Sec_LDA_Class_TrainingSim' num2str(simNum) '_TestingSim' num2str(simulationChoice)],'jpg');
    saveas(appliedF,[save_directory num2str(networkName) SelSims.Name '_Sec_LDA_Class_TrainingSim' num2str(simNum) '_TestingSim' num2str(simulationChoice)],'eps');
end
end
function save_LDA(LDA_Analysis,network,network_load,sim_loaded,currentPath)
cd(currentPath)
save_directory='..\Data\LDA Analysis (Mac)\';
if strcmp(network_load,'z')%Zdenka Code:
    save([save_directory 'Zdenka_' num2str(network.number_of_wires) 'nw_LDA_Analysis_' date],'LDA_Analysis');
elseif strcmp(network_load,'a') %adrian code
    network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
    save([save_directory 'Adrian_' num2str(network.Name) 'LDA_Analysis_' date],'LDA_Analysis','sim_loaded');
end
end

%Graph Functions
function [Graph, binarise_network]=graph_analysis(network,network_load,currentSim,IndexTime,loop)
%% Mac's Analysis: (Graph)

if isempty(IndexTime)
    IndexTime=input(['What Timestamp do you want to analyse? 1-' num2str(size(currentSim.Data,1)) '\n']); %CHOOSE TIMESTAMP
end
%this gives a resistance matrix for the network used for a chosen simulation at a specific timestamp
if loop
    binarise_network='n';
else
    binarise_network=input('Do you want to view the network thresholded/binarised ONLY with low Resistance? \n','s');
end
if binarise_network=='y' %Binarise so we can use Resistance for graph theory analysis
    
    %what we are doing here is creating a matrix of 0 and 1 (high and
    %low resistence), so that we can plot ONLY those nodes/edges that
    %have low resistence and are relevant.
    fprintf('Binarising... \n');
   a= abs(full(currentSim.Data.Rmat{IndexTime}));
    a(isnan(a))=0;
    if network_load=='a'
        
        a(a>0 & a<5000000)=1; %if resistance is low, we make it 1 (on) - need to code it this way because for some reason matlab was showing more than one unique value for each of ON & OFF
        a(a>5000)=0;  %if it is high we make it 0 (off)
    else
        a(a>0 & a<1e7)=1; %if resistance is low, we make it 1 (on)
        a(a>1e4)=0;  %if it is high we make it 0 (off)
    end
    net_mat=a;
    Graph.binarised='Yes - Using Resistance';
    fprintf('Binarisation Complete \n');
else
    fprintf('Binarising... \n');
    net_mat=currentSim.SelLayout.AdjMat; %use standard adjacency matrix
    Graph.binarised='No';
    fprintf('Binarisation Complete \n');
end
%    net_mat=full(simulations(simNum).Data.AdjMat{IndexTime}); %this gives an adj matrix for the network used for a chosen simulation at a specific timestamp

%Global efficiency --> 1/characteristic path length, averaged over the whole network. An estimate of how integrated the network is.
% & Distance Matrix
%The distance matrix contains lengths of shortest paths between all pairs of nodes
[Graph.Distance, Graph.GE] = efficiency_bin(net_mat,0);
fprintf('Global Efficiency Complete \n');

%Topological features will allow us to understand whether it was good or
%bad for classification.

%Local efficiency --> 1/characteristic path length, at each node
[localDistance, Graph.LE]= efficiency_bin(net_mat,1);
fprintf('Local Efficiency Complete \n');

% % Diameter
% Graph.Diameter=diameter(net_mat);


% Density
%Network Density:
% This is defined, for a given set of nodes, as the number of actual edges
% divided by the number of potential edges.
% I.e., the percentage of possible connections that actually exist.
Graph.Density=density_und(net_mat);
fprintf('Density Complete \n');

%Betweeness Centrality
Graph.BC=betweenness_bin(net_mat);
% Node betweenness centrality is the fraction of all shortest paths in the network that contain a given node. Nodes with high values of betweenness centrality participate in a large number of shortest paths.
% fprintf('Betweenness Centrality Complete \n');

%Eigenvector Centrality
Graph.EC=eigenvector_centrality_und(full(net_mat));
% Eigenector centrality is a self-referential measure of centrality -- nodes have high eigenvector centrality if they connect to other nodes that have high eigenvector centrality.
% fprintf('Eigenvector Centrality Complete \n');

%Subgraph Centrality
% Graph.SC=subgraph_centrality(full(net_mat));
%The subgraph centrality of a node is a weighted sum of closed walks of different lengths in the network starting and ending at the node.
% fprintf('Subgraph Centrality Complete \n');

%communicability --> an estimate of the ease with which each pair of nodes can connect in the network
Graph.COMM = expm(net_mat);
fprintf('Communicability Complete \n');

%Clustering Coefficient
[Graph.GlobalClust,Graph.AvgLocalClust, Graph.Clust] = clustCoeff(net_mat);
% Graph.Clust2=clustering_coef_bu(net_mat);
% Graph.GlobalC1ust = number of triangle loops / number of connected triples
% Graph.AvgLocalClust = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)
% Graph.Clust = a 1xN vector of clustering coefficients per node (where mean(C) = C2)
fprintf('Clustering Coefficient Complete \n');


% Graph.Distance=distance_bin(net_mat);
% fprintf('Distance Matrix Complete \n');

%Path Length
Graph.Path = path_length(net_mat);
[Graph.CharPath,Graph.Efficiency,Graph.ecc,Graph.Radius,Graph.Diameter]=charpath(Graph.Distance);
%Average Path Length
Graph.AvgPath=mean(Graph.Path);
fprintf('Path Length Complete \n');
fprintf('Diameter Complete \n');

%Connectivity Magnitude Measure of the network:
%Average number of connections divided by maximum number of possible connections
% In the human brain it is a ratio of ~10^-17 neural connections to possible connections - Kozachkov et al. (2019). How neural circuits achieve and use stable dynamics
Graph.ConnectivityMagnitude=Graph.AvgPath/(length(net_mat)*length(net_mat));

%Small World Propensity:
Graph.SmallWorldProp=small_world_propensity(net_mat);
fprintf('Small World Propensity Complete \n');

%modularity --> an estimate of how segregated the network is
[Graph.Ci,Graph.Q] = community_louvain(net_mat,1);
%The Ci term is the module assignment for each node
%The Q term is the 'quality' of the partition --> how modular the network is.
% -- this should tell us how well we can classify
fprintf('Ci & Q Complete \n');

%Modularity:
Graph.Modularity=modularity_und(full(net_mat));
%The optimal community structure is a subdivision of the network into nonoverlapping groups of nodes in a way that maximizes the number of within-group edges, and minimizes the number of between-group edges.
%The modularity is a statistic that quantifies the degree to which the network may be subdivided into such clearly delineated groups
fprintf('Modularity Complete \n');

%degree --> a count of how many edges are connected to each node
Graph.DEG = degrees_und(net_mat);
fprintf('Degree Complete \n');

%participation coefficient --> an estimate of how integrative a node is
Graph.P = participation_coef(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'
fprintf('Participation Coefficient Complete \n');

%module degree z-score --> an estimate of how segregated a node is
Graph.MZ = module_degree_zscore(net_mat,Graph.Ci);
%Ci from 'community_louvain.m'
fprintf('Module Z Score Complete \n');

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

