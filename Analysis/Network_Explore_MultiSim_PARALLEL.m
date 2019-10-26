%---------------------------
% 22/04/19
% Network Explore Code v1.1
%
% Author: Alon Loeffler
%
%
% THE NEXT SCRIPT TO RUN AFTER THIS IS MultiSims_Graph_Theory.m
%
% IMPORTANT:
% Threshold = Degree of greater than 1
% Binarise Threshold = Only low resistence junctions + wires
% Full Graph = all degrees, all resistences
% --------------------------
% dbstop if error
%% Load Data

function [Explore, threshold, network, SelSims,analysis_type]=Network_Explore_MultiSim_PARALLEL(numNW,simNum,biasType,SelSims)
set(0,'DefaultFigureVisible','off')


computer=getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        currentPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis';
    case '' %if on linux
        currentPath='/headnode2/aloe8475/CODE/ZdenkaCode/MultiElectrodes';
        %             currentPath='/import/silo2/aloe8475/Documents/CODE/Zdenka Code/MultiElectrodes/';
    case 'LAPTOP-S1BV3HR7'
        currentPath='D:\alon_\Research\PhD\CODE\Analysis';
        %case '' %--- Add other computer paths (e.g. Mike)
end
cd(currentPath);
load_data_question='d';
%     clearvars -except load_data_question currentPath simNum numNW
%     close all
%load network data
%     [network, network_load, simulations,sim_loaded,numNetworks, explore_network]= load_data(currentPath);
simDetails=[numNW; simNum];
network_load='z';
% [network,sim_loaded, explore_network, numNetworks]=Load_Zdenka_Code(simDetails,biasType,SelSims);
% dataPath='/import/silo2/aloe8475/Documents/CODE/Data/Raw/Networks/Zdenka Networks/';
computer     = getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        dataPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Zdenka Networks\';
    case ''
        dataPath='/headnode2/aloe8475/CODE/Data/Raw/Networks/Zdenka Networks/';
end
load([dataPath 'AdriantoZdenka' num2str(numNW) 'nw_simulation1.mat'],'SelNet');
network=SelNet;
network.Simulations=[]; %clear previous simulations saved with the network file
network.Simulations=[SelSims];
simulations=network.Simulations;
cd(currentPath);

%% Choose Simulations
%-------
%Choose Network and Simulation for training or exploring

for currentSimulation=1:length(simulations)
    fprintf(num2str(currentSimulation));
    if load_data_question~='a'
        
        networkNum=1;%input(['Which Network # do you want to explore? 1 - ' num2str(length(network)) '\n']);
        %             simNum=currentSimulation;%input(['Which Simulation # do you want to explore? 1 - '  num2str(length(network(networkNum).Simulations)) '\n']); %% CHANGE WHICH SIMULATION YOU WANT TO TEST HERE.
    end
    
    %     if size(simulations,1)==1
    %         currentSim=simulations{simNum};
    %     else
    if strcmp(biasType{1},'TimeDelay')
        currentSim=simulations{currentSimulation};
    else
        currentSim=simulations;
    end
    %     end
    
    % fprintf(['Simulation: ' network(networkNum).Name currentSim.Name ' selected \n\n\n']);
    
    
    %% Analysis:
    i = 1;
    while i == 1
        %Choose Analysis to perform
        %% CHANGE BETWEEN T AND E HERE
        if exist('biasType','var') %if bias type exists
            if length(biasType)>1 %if it's greater than one, we are doing the time-delay LDA
                analysis_type='t';
            else
                if strcmp(biasType,'DCandWait')
                    analysis_type='e';
                elseif strcmp(biasType,'TimeDelay')
                    analysis_type='t'; %lower(input('Which analysis would you like to perform? T = Time b/w Pulses, E = Early/Mid/Late/Never \n','s'));
                end
            end
        end
        if analysis_type=='e'
            %% Exploratary analysis of simulation
            
            if network_load == 'z' %making sure zdenka code is loaded
                
                %Choose a time to Explore Simulation:
                
                %Find pulse centres
                if min(currentSim.Data.VSource1) < max(currentSim.Data.VSource1)
                    V1 = currentSim.Data.VSource1 - min(currentSim.Data.VSource1(currentSim.Data.VSource1>0)); %if AmpOff is greater than 0 (which it should be)
                else
                    V1=currentSim.Data.VSource1;
                end
                pulseEnds    = [];
                pulseStarts  = [];
                
                j = 1;
                
                isPulse = false;
                
                for m = 1:numel(V1)
                    if V1(m) > 0 && ~isPulse
                        pulseStarts(j) = m;
                        isPulse = true;
                    end
                    
                    if V1(m) <= 0 && isPulse
                        pulseEnds(j) = m - 1;
                        isPulse = false;
                        j = j + 1;
                    end
                end
                
                if numel(pulseStarts) > numel(pulseEnds)
                    pulseEnds(j) = numel(V1);
                end
                
                pulseCentres = floor((pulseStarts + pulseEnds)/2);
                if length(pulseCentres)==currentSim.Settings.SetFreq
                    pulseCentres=[pulseCentres pulseEnds(end-1)]; %take the second last pulse time
                end
            end
            %Choose a time to Explore Simulation:
            %Find the centres of when the voltage is on
            if min(currentSim.Data.VSource1) < max(currentSim.Data.VSource1)
                V1 = currentSim.Data.VSource1 - min(currentSim.Data.VSource1);
            else
                V1=currentSim.Data.VSource1;
            end
            pulseEnds    = [];
            pulseStarts  = [];
            
            j = 1;
            
            isPulse = false;
            
            for i = 1:numel(V1)
                if V1(i) > 0 && ~isPulse
                    pulseStarts(j) = i;
                    isPulse = true;
                end
                
                if V1(i) <= 0 && isPulse
                    pulseEnds(j) = i - 1;
                    isPulse = false;
                    j = j + 1;
                end
            end
            
            if numel(pulseStarts) > numel(pulseEnds)
                pulseEnds(j) = numel(V1);
            end
            
            pulseCentres = floor((pulseStarts + pulseEnds)/2);
            if length(pulseCentres)==currentSim.Settings.SetFreq
                pulseCentres=[pulseCentres pulseEnds(end-1)]; %take the second last pulse time
            end
            
            if isempty(pulseCentres)
                fprintf(num2str(currentSimulation));
            end
            MAX_PULSE_CENTRES=11;
            shg
            for time=1:length(pulseCentres)
                if ~isempty(currentSim.Data.Rmat{pulseCentres(time)})
                    [TimeData(time).Explore{currentSimulation},TimeData(time).threshold{currentSimulation}]=explore_simulation(currentSim,network,network_load,simNum,currentPath,currentSimulation,simulations,pulseCentres,time);
                    progressBar(time,length(pulseCentres))
                    shg
                else
                    TimeData(time).Explore{currentSimulation}=[];
                    TimeData(time).threshold{currentSimulation}=[];
                    progressBar(time,length(pulseCentres))
                end
            end
            if length(pulseCentres<11)
                for time=length(pulseCentres)+1:MAX_PULSE_CENTRES
                    TimeData(time).Explore{currentSimulation}=[];
                    TimeData(time).threshold{currentSimulation}=[];
                end
            end
            %% Saving Explore
            if currentSimulation==length(simulations)
                save_state='n';%lower(input('Would you like to save the Exploration Analysis? y or n \n','s'));
                Explore={TimeData.Explore};
                threshold={TimeData.threshold};
                if save_state=='y'
                    Explore={TimeData.Explore};
                    threshold={TimeData.threshold};
                    save_explore(Explore,network(networkNum),network_load,currentPath,simNum,threshold,simulations,analysis_type);
                    i=i+1; %get out of while loop this loop finishes
                end
            end
            i=i+1;
        elseif analysis_type == 't'
            %% ALON 25/09/19
            % BIN PATH LENGTHS ACROSS PAIRINGS OF ELECTRODES FOR TIME DELAY ANALYSIS
            
            % TIME-BASED ANALYSIS USING ZDENKA CODE:
            if network_load == 'z' %making sure zdenka code is loaded
                
                %Choose a time to Explore Simulation:
                
                %Find pulse centres
                if min(currentSim.Data.VSource1) < max(currentSim.Data.VSource1)
                    V1 = currentSim.Data.VSource1 - min(currentSim.Data.VSource1(currentSim.Data.VSource1>0)); %if AmpOff is greater than 0 (which it should be)
                else
                    V1=currentSim.Data.VSource1;
                end
                
                pulseEnds    = [];
                pulseStarts  = [];
                
                j = 1;
                
                isPulse = false;
                
                for m = 1:numel(V1)
                    if V1(m) > 0 && ~isPulse
                        pulseStarts(j) = m;
                        isPulse = true;
                    end
                    
                    if V1(m) <= 0 && isPulse
                        pulseEnds(j) = m - 1;
                        isPulse = false;
                        j = j + 1;
                    end
                end
                
                if numel(pulseStarts) > numel(pulseEnds)
                    pulseEnds(j) = numel(V1);
                end
                
                pulseCentres = floor((pulseStarts + pulseEnds)/2);
                %Find midpoint between the last pulse centre and the second
                %last pulse centre.
                pulseTimeMid=floor((pulseCentres(end)+pulseCentres(end-1))/2); %find the middle of the third pulse and the 4th pulse (i.e. the time difference)
                pulseTimeEnd=pulseStarts(end)-1; %the time just before the last pulse
                pulseCentres=[pulseCentres(1:end-1) pulseTimeMid pulseTimeEnd pulseCentres(end)];
                for time=1:length(pulseCentres) %Alon to change to var
                    [TimeData(time).Explore{currentSimulation},TimeData(time).threshold{currentSimulation}]=explore_simulation(currentSim,network,network_load,simNum,currentPath,currentSimulation,simulations,pulseCentres,time);
                    progressBar(time,length(pulseCentres));
                    
                end
                %% Saving Explore
                %                 if currentSimulation==length(simulations)
                save_state='n';%lower(input('Would you like to save the Exploration Analysis? y or n \n','s'));
                Explore={TimeData.Explore};
                threshold={TimeData.threshold};
                if save_state=='y'
                    Explore={TimeData.Explore};
                    threshold={TimeData.threshold};
                    save_explore(Explore,network(networkNum),network_load,currentPath,simNum,threshold,simulations,analysis_type);
                    i=i+1; %get out of while loop this loop finishes
                end
                %                 end
                i=i+1;
            end
        end
    end
end

%--------------------------------------------------------------------------
end

%% FUNCTIONS

%Loading Functions
% function [network, network_load, simulations, sim_loaded, numNetworks, explore_network] = load_data(currentPath,simNum)
%
% %% Load Data
% %Ask to load Zdenka or Adrian:
% network_load=lower(input('Which Network do you want to analyse? Z - Zdenka, A - Adrian \n','s'));
%
% if strcmp(network_load,'a')
%     %Get current network - Adrian
%     [network,sim_loaded, explore_network, numNetworks]=Load_Adrian_Code();
%     %unpack simulation data into simulation variable
%     if sim_loaded==1
%         if explore_network=='t' %if we have training and testing simulations
%             tempSim=network.Simulations{2};
%             %         tempSim=num2cell(tempSim);
%             %         network.Simulations(2) = [];
%             %         network.Simulations=[network.Simulations tempSim];
%
%             %number of training + number of testing:
%
%             fprintf(['Your Training Simulations are Simulations 1 - ' num2str(network.numTrainingSims) '\n']);
%             fprintf(['Your Testing Simulations are Simulations ' num2str(network.numTrainingSims +1) ' - ' num2str(network.numTestingSims) '\n']);
%
%             fprintf('\n -------------------------- \nStart Analysis: \n');
%
%         end
%         for i = 1:length(network.Simulations)
%             simulations(i)=network.Simulations(i);
%         end
%     else
%         simulations=network.Simulations;
%     end
%     cd(currentPath);
% elseif strcmp(network_load,'z')
%     % Get network - Zdenka:
%     % D:\alon_\Research\PhD\CODE\Zdenka Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data
%
%     % end
% end
% end

%Exploring Functions:
function [Explore,threshold] = explore_simulation(Sim,network,network_load,simNum,currentPath,currentSimulation,simulations,times,time)

% IMPORTANT:
% Threshold = Degree of greater than 1
% Binarise Threshold = Only low resistence junctions + wires
% Full Graph = all degrees, all resistences

[NodeList.String,NodeList.UserData]=GetNodeList(Sim,network_load);
if network_load=='a'
    NodeList.Value=1:height(Sim.Electrodes);
else
    NodeList.Value=1:length(Sim.Electrodes);
end

%% Timeseries View
%Plot Current
% f=figure;
if network_load=='a'
    drainIndex=find(contains(Sim.Electrodes.Name,'Drain'));
    sourceIndex=find(contains(Sim.Electrodes.Name,'Source'));
else
    drainIndex=find(contains({Sim.Electrodes.Name},'Drain'));
    sourceIndex=find(contains({Sim.Electrodes.Name},'Source'));
    
end
if isempty(drainIndex)
    drain_exist=0;
end
if isempty(sourceIndex)
    source_exist=0;
end
for i = 1:length(sourceIndex)
    if ismember(['ISource' num2str(i)], fieldnames(Sim.Data))
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

%     plot(source./sourceV)
%     title('Source')
%     xlabel('Timestamp (0.01sec)')
%     ylabel('Conductance');


IndexTime=times(time); %choose the last timestamp
%input(['What Timestamp do you want to analyse? 1-' num2str(size(Sim.Data,1)) '\n']); %CHOOSE TIMESTAMP

%% Network View
% Function that plots network view of current and resistance
[ Adj, NumEl, Explore] = network_view(Sim,IndexTime, NodeList,network_load);
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
    threshold_network='g';%lower(input('Do you want to plot the entire Graph or the Thresholded Graph (>1 Degree)? g - entire, t - threshold \n','s'));
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
    [f4, f5, f6, G, Adj, Adj2, Explore, highlightElec, new_electrodes] = graph_view_threshold(Sim,Graph,IndexTime,Explore,G, threshold_network, threshold, drain_exist,network_load,node_indices);
else
    [f4, f5, f6, G, Adj, Explore, highlightElec, new_electrodes] = graph_view(Sim,IndexTime,Explore,G, threshold_network,drain_exist,source_exist,node_indices);
end
%% Graph Theory View
% Function that plots graph theory overlayed on graph view of currents
if threshold_network=='t'
    [f7, f8, f9, f10, f11, f12,f13,f14, Explore,sourceElec, drainElec]= graph_theory_explore_threshold(Sim,G,Adj,Adj2, IndexTime,threshold,threshold_network, Explore, Graph, highlightElec, new_electrodes,node_indices,drain_exist,source_exist,network_load);
    %     fprintf('Graph Theory Complete \n');
    
else
    [f7, f8, f9, f10, f11, f12,f13,f14, Explore, sourceElec, drainElec]= graph_theory_explore(Sim,G,Adj,IndexTime,threshold_network, Explore, Graph, highlightElec, new_electrodes,drain_exist,source_exist);
    %     fprintf('Graph Theory Complete \n');
end



%% Save
%Save Variables
Explore.IndexTime=IndexTime;
Explore.Name=Sim.Settings.Name;
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
% if currentSimulation==length(simulations) %if we are in the final simulation
%     save_explore_plots=lower(input('Would you like to save the plots? y or n \n','s'));
%     cd(currentPath)
%     save_directory='..\Data\Figures\Explore Analysis\';
%     network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
%
%     if save_explore_plots=='y'
%         if threshold_network~='t'
%             % NOTE: 05/06 - for simulations created after this date, we need to
%             % change the file names to Sim.Name, and all this info will be in
%             % there already.
%             %         saveas(f,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Timeseries'],'jpg');
%             %         print(f,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Timeseries.pdf']);
%             saveas(f1,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Conductance_Timeseries' num2str(IndexTime)],'jpg');
%             print(f1,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_Conductance_Timeseries' num2str(IndexTime) '.pdf']);
%             saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
%             print(f2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
%             print(f3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_NetworkView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Currents_Timestamp' num2str(IndexTime)],'jpg');
%             print(f4,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
%             print(f5,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Voltage_Timestamp' num2str(IndexTime)],'jpg');
%             print(f6,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Voltage_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'jpg');
%             print(f7,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'jpg');
%             print(f8,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Module-zScore_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'jpg');
%             print(f9,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Connectivity_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f10,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime)],'jpg');
%             print(f10,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f11,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime)],'jpg');
%             print(f11,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f12,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime)],'jpg');
%             print(f12,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f13,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime)],'jpg');
%             print(f13,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f14,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime)],'jpg');
%             print(f14,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime) '.pdf']);
%
%         elseif threshold_network=='t'
%             saveas(f2,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Currents_Timestamp' num2str(IndexTime)],'jpg');
%             print(f2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f3,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
%             print(f3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_NetworkView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f4,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Currents_Timestamp' num2str(IndexTime)],'jpg');
%             print(f4,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Currents_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f5,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Resistance_Timestamp' num2str(IndexTime)],'jpg');
%             print(f5,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Resistance_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f6,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Voltage_Timestamp' num2str(IndexTime)],'jpg');
%             print(f6,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Voltage_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f7,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime)],'jpg');
%             print(f7,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Participant-Coefficient_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f8,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Module-zScore_Timestamp' num2str(IndexTime)],'jpg');
%             print(f8,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Module-zScore_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f9,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Connectivity_Timestamp' num2str(IndexTime)],'jpg');
%             print(f9,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Connectivity_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f10,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime)],'jpg');
%             print(f10,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_Histograms_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f11,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime)],'jpg');
%             print(f11,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_Distances_From_Source_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f12,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime)],'jpg');
%             print(f12,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Current_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f13,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime)],'jpg');
%             print(f13,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Communicability_Timestamp' num2str(IndexTime) '.pdf']);
%             saveas(f14,[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime)],'jpg');
%             print(f14,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_Explore_THRESHOLD_GraphView_ShortestPath_Overlay_Cluster_Timestamp' num2str(IndexTime) '.pdf']);
%         end
%     end
% else
close all
% end

end

function save_explore(Explore,network,network_load,currentPath,simNum,threshold,Sim,analysis_type)
cd(currentPath);
save_directory='..\Data\Explore Analysis\';
network.Name(regexp(network.Name,'[/:]'))=[]; %remove '/' character because it gives us saving problems
if strcmp(network_load,'a') %adrian code
    if analysis_type=='e'
        save([save_directory 'Adrian_' num2str(network.Name) '_Sim_' num2str(simNum) '_Time+Category_Analysis_Timestamp_' num2str(Explore{1}{simNum}.IndexTime) '_' date],'Explore','threshold','network','Sim','analysis_type','-v7.3');
    else
        save([save_directory 'Adrian_' num2str(network.Name) '_Sim_' num2str(simNum) '_Time+Memory_Analysis_Timestamp_' num2str(Explore{1}{simNum}.IndexTime) '_' date],'Explore','threshold','network','Sim','analysis_type','-v7.3');
    end
else
    if analysis_type=='e'
        save([save_directory 'Zdenka' num2str(network.Name) '_Sim_' num2str(simNum) '_Time+Category_Analysis_Timestamp_' num2str(Explore{1}{simNum}.IndexTime) '_' date],'Explore','threshold','network','Sim','analysis_type','-v7.3');
    else
        save([save_directory 'Zdenka' num2str(network.Name) '_Sim_' num2str(simNum) '_Time+Memory_Analysis_Timestamp_' num2str(Explore{1}{simNum}.IndexTime) '_' date],'Explore','threshold','network','Sim','analysis_type','-v7.3');
    end
end
end

%Graph Functions
function [Graph, binarise_network]=graph_analysis(network,network_load,currentSim,IndexTime)
%% Mac's Analysis: (Graph)
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
    %         fprintf('Binarisation Complete \n');
else
    fprintf('Binarising... \n');
    net_mat=currentSim.SelLayout.AdjMat; %use standard adjacency matrix
    Graph.binarised='No';
    %         fprintf('Binarisation Complete \n');
end

net_mat(isnan(net_mat))=0;
%    net_mat=full(simulations(simNum).Data.AdjMat{IndexTime}); %this gives an adj matrix for the network used for a chosen simulation at a specific timestamp
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

%Assortativity
assortativity_bin(net_mat,0);
%save Adj Matrix
Graph.AdjMat=net_mat;

%save selected time
Graph.IndexTime=IndexTime;
end