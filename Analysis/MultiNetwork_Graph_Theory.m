%% Across Network Exploration:
%-------------------------------------------------------------------------

% This script loads Graph Theory analyses performed in Network_Explore for
% all sizes of networks we have simulated. It then creates random and
% ordered networks of the same size (bootstrapped x100), and compares Graph Theoretical
% concepts between our networks and the random/ordered networks.

% Author: Alon Loeffler
% Date: 19/06/2019
% Version: 1.1

% Changelog:
% 22/07/2019 - Added Guimera and Amaral (2005) colored rectangular
% distribution of participant coefficient and module z score
% 26/06/2019 - Added Participant coefficient, Small World Propensity and
% Betweenness centrality measures.%Initialise Paths
%-------------------------------------------------------------------------

dbstop if error
close all
loadDataMulti=lower(input('Would you like to compare multiple iterations of a network? \n','s'));
loadData=lower(input('Would you like to load cross-network data? \n','s'));
if loadData=='y'
    while 1
        if loadDataMulti=='y'
            currMultiNet=input('Which size Network would you like to load iterations from? 100, 500, 1000 or 2000? \n');
            count=1;
        end
        
        currentLocation=pwd;
        computer=getenv('computername');
        switch computer
            case 'W4PT80T2' %if on desktop at uni - Alon
                if loadData=='y' & loadDataMulti=='n'
                    explore_location='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Continuous DC\';
                    savePath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Ordered+Random Data\';
                    fig_dir='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\';
                elseif loadData=='y' & loadDataMulti=='y'
                    explore_location_origin='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Continuous DC\';
                    explore_location=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\MultiNetwork Analysis\' num2str(currMultiNet) 'nw Alternate NWs\'];
                    savePath=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\MultiNetwork Analysis\'];
                    loadPathRandomOrdered='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Ordered+Random Data\';
                    fig_dir=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\'];
                end
            case '' %if on linux
                explore_location='/headnode2/aloe8475/CODE/Data/Explore Analysis/';
                savePath='/headnode2/aloe8475/CODE/Data/Explore Analysis/Multi-Network Data/';
                %                 fig_dir='/headnode2/aloe8475/CODE/Data/Figures/Explore Analysis/Cross-Network Explore/Graph Theory/';
                
            case 'LAPTOP-S1BV3HR7'
                if loadData=='y' & loadDataMulti=='n'
                    explore_location='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\';
                    savePath='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\Multi-Network Data\';
                    fig_dir='D:\alon_\Research\PhD\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\';
                elseif loadData=='y' & loadDataMulti=='y'
                    explore_location_origin='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\';
                    explore_location=['D:\alon_\Research\PhD\CODE\Data\Explore Analysis\' num2str(currMultiNet) 'nw Alternate NWs\'];
                    savePath=['D:\alon_\Research\PhD\CODE\Data\Explore Analysis\Multi-Network Data\'];
                    fig_dir=['D:\alon_\Research\PhD\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\' num2str(currMultiNet) 'nw Alternate NWs\'];
                    %case '' %--- Add other computer paths (e.g. Mike)
                end
        end
        
        if loadData=='y' & loadDataMulti=='n'
            cd(explore_location)
            %Load three explore analyses:
            fprintf('Loading Baseline Data... \n');
            clear e100 e500 e1000 e2000
            % ADRIAN
            %             e100=load([explore_location 'Adrian_Net_Sx20_NoW100_0325-2019_112338__Sim_1_SourceElectrode_6_DrainElectrode_76_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            %             e500=load([explore_location 'Adrian_Net_Sx20_NoW500_0330-2019_111659__Sim_1_SourceElectrode_18_DrainElectrode_430_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            %             e1000=load([explore_location 'Adrian_Net_Sx20_NoW1000_0606-2019_113353__Sim_1_SourceElectrode_32_DrainElectrode_1000_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            %             e2000=load([explore_location 'Adrian_Net_Sx20_NoW2000_0618-2019_125103__Sim_1_SourceElectrode_158_DrainElectrode_1820_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            % ZDENKA
            %             e100=load([explore_location 'Zdenka_Net_Sx20_NoW109_0930-2019_110329_Length_6_Disp_0_Sim_1_Source_99_Drain_36_Explore_Timestamp_200.mat']);
            %             e500=load([explore_location 'Zdenka_Net_Sx20_NoW500_0330-2019_111659_Length_6_Disp_0_Sim_1_Source_498_Drain_435_Explore_Timestamp_200.mat']);
            %             e1000=load([explore_location 'Zdenka_Net_Sx20_NoW1000_0606-2019_113353_Length_6_Disp_0_Sim_1_Source_997_Drain_408_Explore_Timestamp_200.mat']);
            %             e2000=load([explore_location 'Zdenka_Net_Sx20_NoW2000_0618-2019_125103_Length_6_Disp_0_Sim_1_Source_1999_Drain_935_Explore_Timestamp_200.mat']);
            % CHOOSE WHICH TO LOAD
            waitfor(msgbox('Select the 100nw saved data'));
            [FileName,PathName] = uigetfile('*.mat','Select the 100nw saved data');
            count=1;
            f=fullfile(PathName,FileName);
            e100=load(f);
            waitfor(msgbox('Select the 500nw saved data'));
            [FileName,PathName] = uigetfile('*.mat','Select the 500nw saved data');
            count=1;
            f=fullfile(PathName,FileName);
            e500=load(f);
            waitfor(msgbox('Select the 1000nw saved data'));
            [FileName,PathName] = uigetfile('*.mat','Select the 1000nw saved data');
            count=1;
            f=fullfile(PathName,FileName);
            e1000=load(f);
            waitfor(msgbox('Select the 2000nw saved data'));
            [FileName,PathName] = uigetfile('*.mat','Select the 2000nw saved data');
            count=1;
            f=fullfile(PathName,FileName);
            e2000=load(f);
            
            
            
            fprintf('Data Loaded \n');
            break;
        elseif loadData=='y' & loadDataMulti=='y'
            cd(explore_location_origin)
            %Load three explore analyses:
            fprintf('Loading Baseline Data... \n');
            clear e100 e500 e1000 e2000
            %             e100=load([explore_location_origin 'Adrian_Net_Sx20_NoW100_0325-2019_112338__Sim_1_SourceElectrode_6_DrainElectrode_76_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            %             e500=load([explore_location_origin 'Adrian_Net_Sx20_NoW500_0330-2019_111659__Sim_1_SourceElectrode_18_DrainElectrode_430_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            %             e1000=load([explore_location_origin 'Adrian_Net_Sx20_NoW1000_0606-2019_113353__Sim_1_SourceElectrode_32_DrainElectrode_1000_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            %             e2000=load([explore_location_origin 'Adrian_Net_Sx20_NoW2000_0618-2019_125103__Sim_1_SourceElectrode_158_DrainElectrode_1820_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
            e100=load([explore_location_origin 'Zdenka_Net_Sx20_NoW100_0930-2019_110329_Length_6_Disp_0_Sim_1_Source_99_Drain_36_Explore_Timestamp_200.mat']);
            e500=load([explore_location_origin 'Zdenka_Net_Sx20_NoW500_0330-2019_111659_Length_6_Disp_0_Sim_1_Source_498_Drain_435_Explore_Timestamp_200.mat']);
            e1000=load([explore_location_origin 'Zdenka_Net_Sx20_NoW1000_0606-2019_113353_Length_6_Disp_0_Sim_1_Source_997_Drain_408_Explore_Timestamp_200.mat']);
            e2000=load([explore_location_origin 'Zdenka_Net_Sx20_NoW2000_0618-2019_125103_Length_6_Disp_0_Sim_1_Source_1999_Drain_935_Explore_Timestamp_200.mat']);
            
            fprintf('Data Loaded \n');
            
            cd(explore_location);
            fprintf(['Loading Multi Network Data from Network size ' num2str(currMultiNet) '\n']);
            dinfo = dir('Zdenka*.mat');
            for K = 1 : length(dinfo) %Load all files in selected folder
                thisfile = dinfo(K).name;
                destfile = fullfile(explore_location, thisfile);
                switch currMultiNet
                    case 100
                        e100Multi(K) = load(thisfile);
                    case 500
                        e500Multi(K) = load(thisfile);
                    case 1000
                        e1000Multi(K) = load(thisfile);
                    case 2000
                        e2000Multi(K) = load(thisfile);
                end
            end
        end
        loadDataMulti=lower(input('Would you like to load multiple iterations of another sized Network? \n','s'));
        if loadDataMulti=='n'
            loadDataMulti='y';
            break;
        end
    end
end
%% Run Functions:
if loadDataMulti=='n' %if we are doing our normal analysis without multiple networks:
    loop=0;
    ANN=AI(e500,savePath,loop);
    cElegans=cElegansFun();
    human=humanFun();
    AgNW=AgNWFun(e100, e500, e1000, e2000);
    if loadData == 'y'
        [sizeNetwork,randomOrdered,randomOrdered100]=randomOrderedFun(loadPathRandomOrdered,currentLocation,AgNW,loadDataMulti);%e100,e500,e1000,e2000);
    end
    
    %Plot graphs
    plotAll(sizeNetwork,randomOrdered,randomOrdered100, human, e100, e500, e1000, e2000, AgNW,cElegans,ANN,fig_dir)
    
else % if we are looping through
    loop = 1;
    ANN=AI(e500,savePath,loop);
    cElegans=cElegansFun();
    human=humanFun();
    AgNW=AgNWFunMulti(e100Multi, e500Multi, e1000Multi, e2000Multi);    %function for looping
    %     if loadData == 'y'
    [sizeNetwork,randomOrdered,randomOrdered100]=randomOrderedFun(loadPathRandomOrdered,currentLocation,AgNW,loadDataMulti);
    %     end
    
    plotMulti(sizeNetwork,randomOrdered,randomOrdered100, human, e100Multi, e500Multi, e1000Multi, e2000Multi, AgNW,cElegans,ANN,fig_dir)
end

%% TO DO

%% Example AI Graph Analysis (Recurrent Neural Network)
%     function AI()
%     %Create Sample RNN: - NEED TO TALK TO MAC ABOUT THIS
%     %AI.AdjMat=zeros([500 500]);
%     end

%% FUNCTIONS
% Artificial Neural Networks
function ANN=AI(network,savePath,loop)

cd(savePath)
cd('../../../Analysis/Artificial Neural Network Graph Theory/500-Nodes/');
fprintf('Loading Artificial Neural Network Data \n');
if exist('ANNGraphTheory.mat','file')
    load('ANNGraphTheory.mat');
elseif ~exist('ANNGraphTheory.mat','file')
    load('ANN_Adj.mat'); %load untrained data
    
    B=double(ANN_Adj_Mat_Untrained);
    B_Trained=double(ANN_Adj_Mat_Trained);
    %         s=rng(2);
    %
    %         temp=mnrnd(498,[0.02 0.21 0.25 0.25 0.25  0.02]);
    %         temp2=[];
    %         for i =1:length(temp)
    %             temp2(i)=sum(temp(1:i));
    %         end
    %         temp2 = [1 temp2];
    %         B=zeros(498,498);
    %         C=1;
    %         for i = 1:length(temp)-1
    %             for j = temp2(i):temp2(i+1)
    %                 for k = temp2(i+1):temp2(i+2)
    %                     B(j,k)=1;
    %                 end
    %             end
    %         end
    %
    %         B=B+B';
    %
    
    % UNTRAINED
    [ANN.Ci,ANN.Q] = community_louvain(B,1);
    %Clustering Coefficient
    [ANN.GlobalClust,ANN.AvgLocalClust, ANN.Clust] = clustCoeff(B);
    
    %Participation Coefficient & Module z-Score
    ANN.P = participation_coef(B,ANN.Ci);
    ANN.MZ = module_degree_zscore(B, ANN.Ci);
    % Avg + Std PC + MZ:
    ANN.avgP=mean(ANN.P);
    ANN.stdP=std(ANN.P);
    ANN.avgMZ=mean(ANN.MZ);
    ANN.stdMZ=std(ANN.MZ);
    
    %Communicability:
    ANN.COMM = expm(B);
    %Avg COMM
    ANN.avgCOMM=mean(ANN.COMM);
    ANN.stdCOMM=std(ANN.COMM);
    
    
    %Betweenness Centrality:
    [ANN.BC, ANN.normBC]=betweenness_bin(B);
    %Avg:
    ANN.avgBC=mean(ANN.BC);
    ANN.stdBC=std(ANN.BC);
    
    %ANN Path Length
    ANN.Path = path_length(B);
    ANN.AvgPath=mean(ANN.Path);
    G_up=graph(B,'upper');
    G_low=graph(B,'lower');
    G=graph(B);
    
    ANN.Graph=G;
    %Circuit Rank
    ANN.CircuitRank=numedges(G) - (numnodes(G) - 1);
    ANN.SmallWorldProp=small_world_propensity(B);
    
    %TRAINED
    [ANN.Ci_Trained,ANN.Q_Trained] = community_louvain(B_Trained,1,[],'negative_asym');
    %Clustering Coefficient
    [ANN.GlobalClust_Trained,ANN.AvgLocalClust_Trained, ANN.Clust_Trained] = clustCoeff(B_Trained);
    
    %Participation Coefficient & Module z-Score
    ANN.P_Trained = participation_coef(B_Trained,ANN.Ci_Trained);
    ANN.MZ_Trained = module_degree_zscore(B_Trained, ANN.Ci_Trained);
    % Avg + Std PC + MZ:
    ANN.avgP_Trained=mean(ANN.P_Trained);
    ANN.stdP_Trained=std(ANN.P_Trained);
    ANN.avgMZ_Trained=mean(ANN.MZ_Trained);
    ANN.stdMZ_Trained=std(ANN.MZ_Trained);
    
    %Communicability:
    ANN.COMM_Trained = expm(B_Trained);
    %Avg COMM
    ANN.avgCOMM_Trained=mean(ANN.COMM_Trained);
    ANN.stdCOMM_Trained=std(ANN.COMM_Trained);
    
    
    %         %Betweenness Centrality:
    %         [ANN.BC_Trained, ANN.normBC_Trained]=betweenness_bin(B_Trained);
    %         %Avg:
    %         ANN.avgBC_Trained=mean(ANN.BC_Trained);
    %         ANN.stdBC_Trained=std(ANN.BC_Trained);
    %
    %ANN Path Length
    %% NEED TO MAKE SYMMETRIC MATRIX:
    
    %         A=B_Trained;
    %         B_Trained_Undirected=A+A';
    ANN.Path_Trained = path_length(B_Trained);
    ANN.AvgPath_Trained=mean(ANN.Path_Trained);
    %         G_up_Trained=graph(B_Trained,'upper');
    %         G_low_Trained=graph(B_Trained,'lower');
    G_Trained=graph(B_Trained);
    
    ANN.Graph_Trained=G_Trained;
    %Circuit Rank
    ANN.CircuitRank_Trained=numedges(G_Trained) - (numnodes(G_Trained) - 1);
    ANN.SmallWorldProp_Trained=small_world_propensity(B_Trained);
    
    save('ANNGraphTheory.mat','ANN');
end
end
% C-Elegans:
function elegans=cElegansFun()

%Kaiser M, Hilgetag CC (2006) Non-Optimal Component Placement, but Short Processing Paths, due to Long-Distance Projections in Neural Systems. PLoS Computational Biology 7:e95
% Kötter R (2004) Online retrieval, processing, and visualization of primate connectivity data from the CoCoMac database. Neuroinformatics 2: 127-144.
% Choe Y, McCormick BH, Koh W (2004) Network connectivity analysis on the temporally augmented C. elegans web: A pilot study. Society of Neuroscience Abstracts 30:921.9.
cd('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Organic Networks Connectomes\')


cElegans=load('celegans277neurons.mat');
fprintf('Loading cElegans Data \n');

if exist('cElegansGraphTheory.mat','file')
    load('cElegansGraphTheory.mat','elegans');
else
    
    fprintf('No cElegans Data found, creating cElegans Data \n');
    
    %Elegans Matrix (Taken from above sources):
    A=cElegans.celegans277matrix;
    B=A+A'; %creating undirected matrix from directed matrix.
    [elegans.Ci,elegans.Q] = community_louvain(B,1);
    %Clustering Coefficient
    [elegans.GlobalClust,elegans.AvgLocalClust, elegans.Clust] = clustCoeff(B);
    
    %Participation Coefficient & Module z-Score
    elegans.P = participation_coef(B,elegans.Ci);
    elegans.MZ = module_degree_zscore(B, elegans.Ci);
    % Avg + Std PC + MZ:
    elegans.avgP=mean(elegans.P);
    elegans.stdP=std(elegans.P);
    elegans.avgMZ=mean(elegans.MZ);
    elegans.stdMZ=std(elegans.MZ);
    
    %Communicability:
    elegans.COMM = expm(B);
    %Avg COMM
    elegans.avgCOMM=mean(elegans.COMM);
    elegans.stdCOMM=std(elegans.COMM);
    
    
    %Betweenness Centrality:
    [elegans.BC, elegans.normBC]=betweenness_bin(B);
    %Avg:
    elegans.avgBC=mean(elegans.BC);
    elegans.stdBC=std(elegans.BC);
    
    %Random Path Length
    elegans.Path = path_length(B);
    elegans.AvgPath=mean(elegans.Path);
    G_up=graph(B,'upper');
    G_low=graph(B,'lower');
    G=graph(B);
    
    elegans.Graph=G;
    %Circuit Rank
    elegans.CircuitRank=numedges(G) - (numnodes(G) - 1);
    elegans.SmallWorldProp=small_world_propensity(B);
    
    save('cElegansGraphTheory.mat','elegans');
end
end
% Human Graph Analysis
function human=humanFun()
%Cluster Coeff & Path Length
human.GlobalClust=0.53; %Taken from (Achard et al., 2006)
human.AvgPath=2.49; %Taken from (Achard et al., 2006)
%Participation Coefficient & Module z-Score
human.PLocalHubs=0.2; %High PCoeff - approximates from (Power et al., 2013) %High PCoeff = Hubs / Central areas (Power et al., 2013)
human.PConnectorHubs=0.9; %Low PCoeff
human.AvgP=0.502; %Bertolero, Yeo & D'Esposito - 2015
% LOOK INTO THIS 22/07/19 - Doesn't look right:
human.MZHubs=2.6;%[2.5 10] % within module degree probability less than 0.01 is analogous to having a z-score above 2.5 (Guimera and Amaral, 2005, Joyce et al., 2010)
human.MZNonHubs=median([-2:0.5:2.5]);% (Guimera and Amaral, 2005, Joyce et al., 2010)
human.AvgMZ=0.0001; %Bertolero, Yeo & D'Esposito - 2015


%Small world Prop

%Complexity

%Communicability

%Betweenness Centrality

%Circuit Rank
end
% Random and Ordered Graph Analysis
function [sizeNetwork,randomOrdered, randomOrdered100,network]=randomOrderedFun(savePath,currentLocation,AgNW,loadDataMulti)
cd(currentLocation)
i=1;
while 1
    %     sizeNetwork=input('What Size Random/Ordered Networks would you like to create/load? 100, 500, 1000 or 2000? \n');
    createNewRand=lower(input('Would you like to create new Random and Ordered graphs? (Note this will take 4+ Hours) \n','s'));
    loadPath=savePath;
    % save the data for each network & save the network size in a different
    % variable (sizeNetwork).
    %     switch sizeNetwork
    %         case 100
    %             i=1;% if we accidentally load the same network twice, this will override it
    %             network{i}=e100;
    %             sizeNetwork=100;
    %         case 500
    %             i=2;
    %             network{i}=e500;
    %
    %             sizeNetwork=500;
    %         case 1000
    %             i=3;
    %             network{i}=e1000;
    %
    %             sizeNetwork=1000;
    %         case 2000
    %             i=4;
    %             network{i}=e2000;
    %
    %             sizeNetwork=2000;
    %     end
    for k = 1:length(AgNW.NumWires)
        sizeNetwork(k)=AgNW.NumWires(k);
    end
    if createNewRand=='n'
        if exist([loadPath 'Random_Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat'], 'file') == 2 | exist([loadPath 'Random_Ordered_Graphs_networks.mat'], 'file') == 2 %&& exist([loadPath 'Random_Graphs_' num2str(sizeNetwork) 'nw.mat'],'file') == 2 %2 because .mat file
            
            fprintf(['Loading Ordered and Random Graphs (' num2str(sizeNetwork) 'nodes)... \n \n']);
            %                     load([loadPath 'Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat']);
            %                     load([loadPath 'Random_Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat']);
            load([loadPath 'Random_Ordered_Graphs_networks.mat']);
            break;
            
            
        else
            fprintf('Ordered and Random Graphs have not been created yet \n');
            fprintf('Creating New Graphs \n');
            %             [random, ordered, randomOrdered100, ordered100, Parameters]=createRandom_Ordered_Graphs(sizeNetwork,AgNW);
            if loadDataMulti=='n'
                [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW);
            else
                [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs_Multi(AgNW);
            end
            
        end
    else
        if exist([loadPath 'Random_Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat'], 'file') == 2 | exist([loadPath 'Random_Ordered_Graphs_ALL_networks.mat'], 'file') == 2 %%&& exist([loadPath 'Random_Graphs_' num2str(sizeNetwork) 'nw.mat'],'file') == 2 %2 because .mat file
            overwrite=lower(input(['Graphs (' num2str(sizeNetwork) 'nodes) have already been created, would you like to overwrite? \n'],'s'));
            if overwrite =='n'
                fprintf(['Loading Ordered and Random Graphs (' num2str(sizeNetwork) 'nodes) \n']);
                %                 load([loadPath 'Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat']);
                %                 load([loadPath 'Random_Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat']);
                load([loadPath 'Random_Ordered_Graphs_ALL_networks.mat']);
            else
                fprintf('Creating New Graphs \n');
                %                 [random, ordered, randomOrdered100, ordered100, Parameters]=createRandom_Ordered_Graphs(sizeNetwork,AgNW);
                if loadDataMulti=='n'
                    [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW);
                else
                    [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs_Multi(AgNW);
                end
                fprintf('New Graphs Created \n');
            end
        else
            fprintf('Creating New Graphs \n');
            if loadDataMulti=='n'
                [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW);
            else
                [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs_Multi(AgNW);
            end
            fprintf('New Graphs Created \n');
        end
    end
    %save in struct so we can load multiple networks at once:
    load2=lower(input('Would you like to load another network? \n','s'));
    if load2=='n'
        break;
    end
end

fprintf('Graphs Loaded... \n Creating Plots... \n');

%% Random & Ordered
%Number of Nodes
% for z=1:length(Net)
%     while isempty(sizeNetwork) %if we only load 500, 1000 or 2000, we need to increase our counter to skip over the empty arrays.
%         z=z+1;
%     end
% for i = 1:4
% %     for j =1:11
% %         randomOrdered.numNodes(i,j)=height(randomOrdered.Graph{i,j}.Nodes);
% %         %         ordered(i).numNodes=height(ordered(i).Graph.Nodes);
% %         randomOrdered.numEdges(i,j)=height(randomOrdered.Graph{i,j}.Edges);
% %         %         ordered(i).numEdges=height(ordered(i).Graph.Edges);
% %     end
% %
% %     randomOrdered100.AvgNodes=mean([randomOrdered.numNodes{i,:}]);
% %     %     ordered100.AvgNodes=mean([ordered(:).numNodes]);
% %
% %     %Number of Edges
% %     randomOrdered100.AvgEdges=mean([randomOrdered.numEdges{i,:}]);
%     %     ordered100.AvgEdges=mean([ordered(:).numEdges]);
%
%     % participation coefficient (mean)
%
%     % avg and std pcoeff across 100 bootraps per node:
%     randomOrdered100.AvgPCoeff(i,:)=mean([randomOrdered.P{i,:}],2);
%     %     ordered100.AvgPCoeff=mean([ordered(:).P],2);
%     randomOrdered100.StdPCoeff(i,:)=std([randomOrdered.P{i,:}],[],2);
%     %     ordered100.StdPCoeff=std([ordered(:).P],[],2);
%
%     % Avg & std PCoeff (across 100 bootstraps) and across nodes:
%     randomOrdered100.AvgAvgPCoeff(i,:)=mean(mean([randomOrdered.P{i,:}]),2);
%     %     randomOrdered100.StdAvgPCoeff=std(mean([random(:).P]),[],2);
%     %     ordered100.AvgAvgPCoeff=mean(mean([randomOrdered(:).P]),2);
%     %     ordered100.StdAvgPCoeff=std(mean([ordered(:).P]),[],2);
%     %Module z-Score
%     % avg and std mz across 100 bootraps per node:
%     randomOrdered100.AvgMZ(i,:)=mean([randomOrdered.MZ{i,:}],2);
%     %     ordered100.AvgMZ=mean([ordered(:).MZ],2);
%     randomOrdered100.StdMZ(i,:)=std([randomOrdered.MZ{i,:}],[],2);
%     %     ordered100.StdMZ=std([ordered(:).MZ],[],2);
%
%     %avg and std mz across 100 bootstraps and across nodes:
%     randomOrdered100.AvgAvgMZ(i,:)=mean(mean([randomOrdered.MZ{i,:}]),2);
%     randomOrdered100.StdAvgMZ(i,:)=std(mean([randomOrdered.MZ{i,:}]),[],2);
%     %     ordered100.AvgAvgMZ=mean(mean([ordered(:).MZ]),2);
%     %     ordered100.StdAvgMZ=std(mean([ordered(:).MZ]),[],2);
%     %Small World Prop
%     randomOrdered100.AvgSmallWorldProp(i,:)=mean([randomOrdered.SmallWorldProp{i,:}],2);
%     randomOrdered100.StdSmallWorldProp(i,:)=std([randomOrdered.SmallWorldProp{i,:}],[],2);
%     %     ordered100.AvgSmallWorldProp=mean([ordered(:).SmallWorldProp],2);
%     %     ordered100.StdSmallWorldProp=std([ordered(:).SmallWorldProp],[],2);
%
%     % communicability
%     randomOrdered100.AvgCOMM(i,:)=mean([randomOrdered.COMM{i,:}],2);
%     randomOrdered100.StdCOMM(i,:)=std([randomOrdered.COMM{i,:}],[],2);
%     % complexity (from the Sporns, Tononi and Edelman paper I sent through the other day -- still waiting on code from Olaf, but will send through when I get it).
%
%     % betweeness_centrality
%
%     %Average Normalised Betweenness (BC/[(N-1)(N-2)])
%     %     randomOrdered100.AvgNormBC=mean([random(:).normBC]);
%     %     randomOrdered100.StdNormBC=std([random(:).normBC]);
%     %     ordered100.AvgNormBC=mean([ordered(:).normBC]);
%     %     ordered100.StdNormBC=std([ordered(:).normBC]);
% end
end
% Silver Nanowire Networks Analysis
function AgNW=AgNWFun(e100, e500, e1000, e2000)
%% AgNW
%Circuit Rank
AgNW.CircuitRank=[e100.Explore.GraphTheory.CircuitRank e500.Explore.GraphTheory.CircuitRank e1000.Explore.GraphTheory.CircuitRank e2000.Explore.GraphTheory.CircuitRank];
AgNW.GlobalClust=[e100.Explore.GraphTheory.GlobalClust, e500.Explore.GraphTheory.GlobalClust, e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
%Degree
AgNW.Degree={e100.Explore.GraphTheory.DEG; e500.Explore.GraphTheory.DEG; e1000.Explore.GraphTheory.DEG; e2000.Explore.GraphTheory.DEG};
AgNW.NumWires=[length(e100.Explore.GraphView.AdjMat),length(e500.Explore.GraphView.AdjMat),length(e1000.Explore.GraphView.AdjMat),length(e2000.Explore.GraphView.AdjMat)];
%Path
AgNW.AvgPath=[e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
%Small World Prop
AgNW.SmallWorldProp=[e100.Explore.GraphTheory.SmallWorldProp, e500.Explore.GraphTheory.SmallWorldProp, e1000.Explore.GraphTheory.SmallWorldProp e2000.Explore.GraphTheory.SmallWorldProp];
%Communicability
AgNW.AvgCOMM=[mean(mean(e100.Explore.GraphTheory.COMM)) mean(mean(e500.Explore.GraphTheory.COMM)) mean(mean(e1000.Explore.GraphTheory.COMM)) mean(mean(e2000.Explore.GraphTheory.COMM))];
AgNW.StdCOMM=[std(mean(e100.Explore.GraphTheory.COMM)) std(mean(e500.Explore.GraphTheory.COMM)) std(mean(e1000.Explore.GraphTheory.COMM)) std(mean(e2000.Explore.GraphTheory.COMM))];
%Betweenness Centrality
AgNW.AvgBC=[mean(e100.Explore.GraphTheory.BC) mean(e500.Explore.GraphTheory.BC)  mean(e1000.Explore.GraphTheory.BC) mean(e2000.Explore.GraphTheory.BC)];
AgNW.StdBC=[std(e100.Explore.GraphTheory.BC) std(e500.Explore.GraphTheory.BC)  std(e1000.Explore.GraphTheory.BC) std(e2000.Explore.GraphTheory.BC)];
save('AgNWGraphTheory.mat','AgNW');
end
function AgNW=AgNWFunMulti(e100Multi, e500Multi, e1000Multi, e2000Multi)
%100 Nw

for z = 1:length(e100Multi.Explore)
    GlobalClust(z)=e100Multi.Explore{z}.GraphTheory.GlobalClust;
    
    AvgPath(z)=e100Multi.Explore{z}.GraphTheory.AvgPath;
    CircuitRank(z)=e100Multi.Explore{z}.GraphTheory.CircuitRank;
    SmallWorldProp(z)=e100Multi.Explore{z}.GraphTheory.SmallWorldProp;
    %         BC(z,:)=e100Multi.Explore{z}.GraphTheory.BC;
    PCoeff{z}=e100Multi.Explore{z}.GraphTheory.P;
    MZ{z}=e100Multi.Explore{z}.GraphTheory.MZ;
    %Degree
    Degree{z}=e100Multi.Explore{z}.GraphTheory.DEG;
end
e100.AvgGlobalClust=mean([GlobalClust]);
e100.StdGlobalClust=std([GlobalClust])./sqrt(length([GlobalClust]));
e100.AvgDEG=mean([Degree{:}]);
e100.StdDEG=std([Degree{:}])./sqrt(length([Degree{z}]));
e100.AvgAvgPath=mean([AvgPath]);
e100.StdAvgPath=std([AvgPath])./sqrt(length([AvgPath]));
e100.AvgCircuitRank=mean([CircuitRank]);
e100.StdCircuitRank=std([CircuitRank]);
e100.AvgSmallWorldProp=mean([SmallWorldProp]);
e100.StdSmallWorldProp=std([SmallWorldProp]);
%     e100.AvgBC=mean([BC]);
%     e100.StdBC=std([BC]);
e100.AvgPCoeff=mean(vertcat(PCoeff{:}));
e100.StdPCoeff=std(vertcat(PCoeff{:}));
e100.AvgMZ=mean(vertcat(MZ{:}));
e100.StdMZ=std(vertcat(MZ{:}));
clear PCoeff SmallWorldProp GlobalClust AvgPath CircuitRank MZ Degree%BC
%500 Nw
for z = 1:length(e500Multi.Explore)
    GlobalClust(z)=e500Multi.Explore{z}.GraphTheory.GlobalClust;
    AvgPath(z)=e500Multi.Explore{z}.GraphTheory.AvgPath;
    CircuitRank(z)=e500Multi.Explore{z}.GraphTheory.CircuitRank;
    SmallWorldProp(z)=e500Multi.Explore{z}.GraphTheory.SmallWorldProp;
    %         BC(z,:)=e500Multi.Explore{z}.GraphTheory.BC;
    PCoeff{z}=e500Multi.Explore{z}.GraphTheory.P;
    MZ{z}=e500Multi.Explore{z}.GraphTheory.MZ;
    %Degree
    Degree{z}=e500Multi.Explore{z}.GraphTheory.DEG;
end
e500.AvgGlobalClust=mean([GlobalClust]);
e500.StdGlobalClust=std([GlobalClust])./sqrt(length([GlobalClust]));
e500.AvgDEG=mean([Degree{:}]);
e500.StdDEG=std([Degree{:}])./sqrt(length([Degree{:}]));
e500.AvgAvgPath=mean([AvgPath]);
e500.StdAvgPath=std([AvgPath])./sqrt(length([AvgPath]));
e500.AvgCircuitRank=mean([CircuitRank]);
e500.StdCircuitRank=std([CircuitRank]);
e500.AvgSmallWorldProp=mean([SmallWorldProp]);
e500.StdSmallWorldProp=std([SmallWorldProp]);
%     e500.AvgBC=mean([BC]);
%     e500.StdBC=std([BC]);
e500.AvgPCoeff=mean(vertcat(PCoeff{:}));
e500.StdPCoeff=std(vertcat(PCoeff{:}));
e500.AvgMZ=mean(vertcat(MZ{:}));
e500.StdMZ=std(vertcat(MZ{:}));
clear PCoeff SmallWorldProp GlobalClust AvgPath CircuitRank MZ Degree%BC
%1000 Nw
for z = 1:length(e1000Multi.Explore)
    GlobalClust(z)=e1000Multi.Explore{z}.GraphTheory.GlobalClust;
    AvgPath(z)=e1000Multi.Explore{z}.GraphTheory.AvgPath;
    CircuitRank(z)=e1000Multi.Explore{z}.GraphTheory.CircuitRank;
    SmallWorldProp(z)=e1000Multi.Explore{z}.GraphTheory.SmallWorldProp;
    %         BC(z,:)=e1000Multi.Explore{z}.GraphTheory.BC;
    PCoeff{z}=e1000Multi.Explore{z}.GraphTheory.P;
    MZ{z}=e1000Multi.Explore{z}.GraphTheory.MZ;
    %Degree
    Degree{z}=e1000Multi.Explore{z}.GraphTheory.DEG;
end
e1000.AvgGlobalClust=mean([GlobalClust]);
e1000.StdGlobalClust=std([GlobalClust])./sqrt(length([GlobalClust]));
e1000.AvgDEG=mean([Degree{:}]);
e1000.StdDEG=std([Degree{:}])./sqrt(length([Degree{:}]));
e1000.AvgAvgPath=mean([AvgPath]);
e1000.StdAvgPath=std([AvgPath])./sqrt(length([AvgPath]));
e1000.AvgCircuitRank=mean([CircuitRank]);
e1000.StdCircuitRank=std([CircuitRank]);
e1000.AvgSmallWorldProp=mean([SmallWorldProp]);
e1000.StdSmallWorldProp=std([SmallWorldProp]);
%     e1000.AvgBC=mean([BC]);
%     e1000.StdBC=std([BC]);
e1000.AvgPCoeff=mean(vertcat(PCoeff{:}));
e1000.StdPCoeff=std(vertcat(PCoeff{:}));
e1000.AvgMZ=mean(vertcat(MZ{:}));
e1000.StdMZ=std(vertcat(MZ{:}));
clear PCoeff SmallWorldProp GlobalClust AvgPath CircuitRank MZ Degree%BC
%2000 Nw
for z = 1:length(e2000Multi.Explore)
    GlobalClust(z)=e2000Multi.Explore{z}.GraphTheory.GlobalClust;
    AvgPath(z)=e2000Multi.Explore{z}.GraphTheory.AvgPath;
    CircuitRank(z)=e2000Multi.Explore{z}.GraphTheory.CircuitRank;
    SmallWorldProp(z)=e2000Multi.Explore{z}.GraphTheory.SmallWorldProp;
    %         BC(z,:)=e2000Multi.Explore{z}.GraphTheory.BC;
    PCoeff{z}=e2000Multi.Explore{z}.GraphTheory.P;
    MZ{z}=e2000Multi.Explore{z}.GraphTheory.MZ;
    %Degree
    Degree{z}=e2000Multi.Explore{z}.GraphTheory.DEG;
end
e2000.AvgGlobalClust=mean([GlobalClust]);
e2000.StdGlobalClust=std([GlobalClust])./sqrt(length([GlobalClust]));
e2000.AvgDEG=mean([Degree{:}]);
e2000.StdDEG=std([Degree{:}])./sqrt(length([Degree{:}]));
e2000.AvgAvgPath=mean([AvgPath]);
e2000.StdAvgPath=std([AvgPath])./sqrt(length([AvgPath]));
e2000.AvgCircuitRank=mean([CircuitRank]);
e2000.StdCircuitRank=std([CircuitRank]);
e2000.AvgSmallWorldProp=mean([SmallWorldProp]);
e2000.StdSmallWorldProp=std([SmallWorldProp]);
%     e2000.AvgBC=mean([BC]);
%     e2000.StdBC=std([BC]);
e2000.AvgPCoeff=mean(vertcat(PCoeff{:}));
e2000.StdPCoeff=std(vertcat(PCoeff{:}));
e2000.AvgMZ=mean(vertcat(MZ{:}));
e2000.StdMZ=std(vertcat(MZ{:}));

%Combine all
AgNW.AvgGlobalClust=[e100.AvgGlobalClust e500.AvgGlobalClust e1000.AvgGlobalClust e2000.AvgGlobalClust];
AgNW.StdGlobalClust=[e100.StdGlobalClust e500.StdGlobalClust e1000.StdGlobalClust e2000.StdGlobalClust];
AgNW.AvgDegree=[e100.AvgDEG e500.AvgDEG e1000.AvgDEG e2000.AvgDEG];
AgNW.StdDegree=[e100.StdDEG e500.StdDEG e1000.StdDEG e2000.StdDEG];
AgNW.AvgAvgPath=[e100.AvgAvgPath e500.AvgAvgPath e1000.AvgAvgPath e2000.AvgAvgPath];
AgNW.StdAvgPath=[e100.StdAvgPath e500.StdAvgPath e1000.StdAvgPath e2000.StdAvgPath];
AgNW.AvgCircuitRank=[e100.AvgCircuitRank e500.AvgCircuitRank e1000.AvgCircuitRank e2000.AvgCircuitRank];
AgNW.StdCircuitRank=[e100.StdCircuitRank e500.StdCircuitRank e1000.StdCircuitRank e2000.StdCircuitRank];
AgNW.AvgSmallWorldProp=[e100.AvgSmallWorldProp e500.AvgSmallWorldProp e1000.AvgSmallWorldProp e2000.AvgSmallWorldProp];
AgNW.StdSmallWorldProp=[e100.StdSmallWorldProp e500.StdSmallWorldProp e1000.StdSmallWorldProp e2000.StdSmallWorldProp];
%     AgNW.AvgBC=[e100.AvgBC e500.AvgBC e1000.AvgBC e2000.AvgBC];
%     AgNW.StdBC=[e100.StdBC e500.StdBC e1000.StdBC e2000.StdBC];
AgNW.AvgPCoeff=[e100.AvgPCoeff e500.AvgPCoeff e1000.AvgPCoeff e2000.AvgPCoeff];
AgNW.StdPCoeff=[e100.StdPCoeff e500.StdPCoeff e1000.StdPCoeff e2000.StdPCoeff];
AgNW.AvgMZ=[e100.AvgMZ e500.AvgMZ e1000.AvgMZ e2000.AvgMZ];
AgNW.StdMZ=[e100.StdMZ e500.StdMZ e1000.StdMZ e2000.StdMZ];
%NEED TO CHANGE LENGTH OF WIRES HERE:
AgNW.NumWires=[100 500 1000 2000];%[length(e100Multi.Explore{1}.GraphView.AdjMat) length(e500Multi.Explore{1}.GraphView.AdjMat) length(e1000Multi.Explore{1}.GraphView.AdjMat) length(e2000Multi.Explore{1}.GraphView.AdjMat)];
end

%% Plot:
function plotAll(Net, randomOrdered100, human, e100, e500, e1000, e2000, AgNW,cElegans,ANN,fig_dir)

% DEGREE
fDeg=figure;
for i = 1:length(AgNW.Degree)
    subplot(2,2,i);
    h=histogram(AgNW.Degree{i},11);
    title([num2str(AgNW.NumWires(i)) 'nw']);
    xlabel('Degree');
    ylabel('Frequency');
end
sgtitle('Degree Distribution Comparison of Nanowire Networks')

%% WATTS STROGATZ
% Small World Analysis

% loadPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Watson Strogatz\';
% b=load([loadPath 'beta.mat']);
% watStr.b=b.beta;
% load([loadPath 'cc.mat']);
% watStr.cc=cc;
% load([loadPath 'pl.mat']);
% watStr.pl=pl;

fWatts=figure('Position',[0 0 1920 1080]);

%Watss-Strogatz Colors:
lightblue=rgb('cyan');
blue=rgb('deep blue');
lightred=rgb('rose pink');
red=rgb('burnt red');


rdif=[blue(1)-lightblue(1)];
r2=blue(1):(-rdif)/100:lightblue(1);
gdif=[blue(2)-lightblue(2)];
g2=blue(2):(-gdif)/100:lightblue(2);
bdif=[blue(3)-lightblue(3)];
b2=blue(3):(-bdif)/100:lightblue(3);

rdif=[red(1)-lightred(1)];
r3=red(1):(-rdif)/100:lightred(1);
gdif=[red(2)-lightred(2)];
g3=red(2):(-gdif)/100:lightred(2);
bdif=[red(3)-lightred(3)];
b3=red(3):(-bdif)/100:lightred(3);

clrs2={[r2; g2; b2]'};
clrs3={[r3; g3; b3]'};

% s(1)=scatter(watStr.pl,mean(watStr.cc),[],[r3;g3;b3;]');
hold on

% f=figure;
while 1
    plotNet=input('Which Network size would you like to plot? 100, 500, 1000, or 2000? \n');
    switch plotNet
        case 100
            plotNet=1;
        case 500
            plotNet=2;
        case 1000
            plotNet=3;
        case 2000
            plotNet=4;
    end
    %     if plotNet<3 %if we only want 100 & 500 networks
    y=[randomOrdered100.AvgGlobalClust{plotNet,:} ANN.GlobalClust ANN.GlobalClust_Trained cElegans.GlobalClust human.GlobalClust  e100.Explore.GraphTheory.GlobalClust e500.Explore.GraphTheory.GlobalClust e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
    %ordered100.AvgGlobalClust
    %ordered100.AvgPath
    x=[randomOrdered100.AvgPath{plotNet,:} ANN.AvgPath ANN.AvgPath_Trained cElegans.AvgPath human.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
    %                p=gscatter(x,y);
    %     hold on
    for i=1:length(x)
        %           s(i+1)=scatter(nanmean(categories.Clust{i}),categories.PathLength{i,j},[],clrs2{1}(i,:));
        s(i)=scatter(x(:,i),y(:,i),[],clrs2{1}(i,:));
        hold on
        e=errorbar(x(1), y(1),randomOrdered100.StdPath);
        e2=errorbar(x(1), y(1),randomOrdered100.StdGlobalClust,'horizontal');
        %             e3=errorbar(x(5), y(5),ordered100.StdPath);
        %             e4=errorbar(x(5), y(5),ordered100.StdGlobalClust,'horizontal');
        e.Color=clrs2{1}(i,:);
        e2.Color=clrs2{1}(i,:);
        %             e3.Color=clrs2{1}(i,:);
        %             e4.Color=clrs2{1}(i,:);
    end
    text(x,y,{[num2str(sizeNetwork) 'node Watts-Strogats'],'500node ANN Untrained','500node ANN Trained','C. Elegans Nw', 'Human Nw (mm-scale regions)','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
    
    %     hold on
    %     e=errorbar(x(1), y(1),randomOrdered100.StdPath);
    %     e2=errorbar(x(1), y(1),randomOrdered100.StdGlobalClust,'horizontal');
    
    %[num2str(sizeNetwork) 'node Ordered Nw']
    
    %     else %if we want all 100, 500, 1000 and 2000 networks plotted
    %         y=[ ANN.GlobalClust ANN.GlobalClust_Trained cElegans.GlobalClust human.GlobalClust  e100.Explore.GraphTheory.GlobalClust e500.Explore.GraphTheory.GlobalClust e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
    %         x=[ANN.AvgPath ANN.AvgPath_Trained cElegans.AvgPath human.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
    %         % xlim([0.05 0.6])
    %         % ylim([2 16])
    %         p=gscatter(x,y);
    %
    %         hold on
    %         text(x,y,{'500node ANN Untrained','500node ANN Trained','C. Elegans Nw', 'Human Nw (mm-scale regions)','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
    %     end
    ylabel('Global Clustering Coefficient');
    xlabel('Global Mean Path Length');
    p.MarkerEdgeColor='b';
    p(:,1).MarkerEdgeColor='r';
    p.LineWidth=1.5;
    
    plot2=lower(input('Would you like to plot another network size? \n','s'));
    if plot2=='n'
        break;
    end
end

fprintf('Figure 1 Complete \n');

%% Small World Log
% f1=figure;
% logx=log10(x);
% logy=log10(y);
% p1=gscatter(logx,y);
% hold on
% hold on
% e=errorbar(logx(1), y(1),randomOrdered100.StdPath);
% e2=errorbar(logx(1), y(1),randomOrdered100.StdGlobalClust);
% errorbar(logx(3), y(3),ordered100.StdPath);
% errorbar(logx(3), y(3),ordered100.StdGlobalClust);
% % xlim([0.05 0.6])
% % ylim([2 16])
% %     ylim([0,max(y)]);
% %     xlim([0,max(logx)]);
% text(logx,y,{'500node Random Nw','Human Nw','500node Ordered Nw','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
% ylabel('Global Clustering Coefficient');
% xlabel('Log10 Global Mean Path Length');
% p1.MarkerEdgeColor='b';
% p1(:,1).MarkerEdgeColor='r';
% p1.LineWidth=1.5;

smallWorldPropRandom=[];
smallWorldPropOrdered=[];

for i = 1:plotNet
    %Circuit Rank:
    smallWorldPropOrderedRandom=[smallWorldPropOrderedRandom randomOrdered100.AvgSmallWorldProp];
    %     smallWorldPropOrdered=[smallWorldPropOrdered ordered100.AvgSmallWorldProp];
    randomLabel{i}=[num2str(sizeNetwork) ' Watts-Strogatz Nw'];
    %     orderedLabel{i}=[num2str(sizeNetwork) ' Ordered Nw'];
end

%Small World Prop
f2=figure;
%smallWorldPropOrdered
x=[smallWorldPropOrderedRandom  ANN.SmallWorldProp cElegans.SmallWorldProp AgNW.SmallWorldProp];
p2=bar(x);
hold on
for i = 1:plotNet
    e=errorbar(x(1), randomOrdered100.StdSmallWorldProp);
    %     e2=errorbar(x(2),ordered100.StdSmallWorldProp);
    % xlim([0.05 0.6])
    % ylim([2 16])
    hold on
end
%orderedLabel
xticklabels([randomLabel,'500node ANN','C. Elegans Nw','100nw','500nw','1000nw','2000nw']);

set(gca, 'XTickLabelRotation', 45)

ylabel('Small World Prop');


fprintf('Figure 2 Complete \n');
%%

f3=figure;
circuitRankRandom=[];
circuitRankOrdered=[];
for i = 1:plotNet
    %Circuit Rank:
    circuitRankRandom=[circuitRankRandom randomOrdered(1).CircuitRank];
    %     circuitRankOrdered=[circuitRankOrdered ordered(1).CircuitRank];
end

circuitRank=[circuitRankRandom circuitRankOrdered ANN.CircuitRank AgNW.CircuitRank];
p3=bar(circuitRank);
% orderedLabel
xticklabels([randomLabel  {'500node ANN','100nw', '500nw', '1000nw','2000nw'}]);
set(gca, 'XTickLabelRotation', 45)

ylabel('Circuit Rank');
hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles
for i=1:length(p3)  % iterate over number of bar objects
    hT=[hT text(p3(i).XData+p3(i).XOffset,p3(i).YData,num2str(p3(i).YData.'), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end
hold on
%Average Participation Coefficient & Module z-Score
fprintf('Figure 3 Complete \n');
%%
f4=figure;

PRandom=[];
POrdered=[];
MZRandom=[];
MZOrdered=[];
stdPRandom=[];
stdPOrdered=[];
stdMZRandom=[];
stdMZOrdered=[];

for i = 1:plotNet
    %Circuit Rank:
    PRandom=[PRandom randomOrdered100.AvgAvgPCoeff];
    %     POrdered=[POrdered ordered100.AvgAvgPCoeff];
    stdPRandom=[stdPRandom randomOrdered100.StdAvgPCoeff];
    %     stdPOrdered=[stdPOrdered ordered100.StdAvgPCoeff];
    MZRandom=[MZRandom randomOrdered100.AvgAvgMZ];
    %     MZOrdered=[MZOrdered ordered100.AvgAvgMZ];
    stdMZRandom=[stdMZRandom randomOrdered100.StdAvgMZ];
    %     stdMZOrdered=[stdMZOrdered ordered100.StdAvgMZ];
end
% POrdered stdMZOrdered MZOrdered

PCoeff=[PRandom  ANN.avgP cElegans.avgP human.AvgP human.PLocalHubs human.PLocalHubs human.PConnectorHubs human.PConnectorHubs  mean(e100.Explore.GraphTheory.P), mean(e500.Explore.GraphTheory.P), mean(e1000.Explore.GraphTheory.P) mean(e2000.Explore.GraphTheory.P)];
stdPCoeff=[stdPRandom  ANN.avgP cElegans.avgP, 0, 0, 0, 0, 0,  std(e100.Explore.GraphTheory.P), std(e500.Explore.GraphTheory.P), std(e1000.Explore.GraphTheory.P) std(e2000.Explore.GraphTheory.P)];
stdMZ=[stdMZRandom ANN.stdMZ cElegans.stdMZ,0,0,0,0,0, std(e100.Explore.GraphTheory.MZ), std(e500.Explore.GraphTheory.MZ), std(e1000.Explore.GraphTheory.MZ) std(e2000.Explore.GraphTheory.MZ)];
MZ=[MZRandom  ANN.avgMZ cElegans.avgMZ human.AvgMZ human.MZHubs human.MZNonHubs human.MZHubs human.MZNonHubs  mean(e100.Explore.GraphTheory.MZ), mean(e500.Explore.GraphTheory.MZ), mean(e1000.Explore.GraphTheory.MZ) mean(e2000.Explore.GraphTheory.MZ)];

p4=gscatter(PCoeff,MZ);
% e=errorbar(MZ, stdMZ);
% e2=errorbar(PCoeff,stdPCoeff);

%High PCoeff = Hubs / Central areas (Power et al., 2013)
%stdMZRandom
text(PCoeff,MZ,[randomLabel, '500node ANN', 'C. Elegans Nw', 'Human Average', 'Human Connector Local Provincial Hub','Human Local Peripheral Node','Human Connector Hub','Human Satellite Connector', '100nw Avg', '500nw Avg', '1000nw Avgs','2000nw Avg'],'NorthWest');
xlabel('Average Participant Coefficient Coefficient');
ylabel('Average Module z-Score');
p4.MarkerEdgeColor='b';
p4(:,1).MarkerEdgeColor='r';
p4.LineWidth=1.5;
hold on

fprintf('Figure 4 Complete \n');

%% Plot Guimera & Amaral rectangles:
f5=figure;
for i=1:plotNet
    subplot(2,2,i)
    if i == 1
        network{i}=e100;
    elseif i==2
        network{i}=e500;
    elseif i==3
        network{i}=e1000;
    elseif i==4
        network{i}=e2000;
    end
    guimera(network{i},randomOrdered100,i); %change network here
    fprintf(['Figure 5 part ' num2str(i) ' Complete \n']);
end

fprintf('Figure 5 Complete \n');

% %Communicability
% f6=figure;
% COMM=[randomOrdered100.AvgCOMM(1) ordered100.AvgCOMM(1) AgNW.AvgCOMM];
% stdCOMM=[0 0 AgNW.StdCOMM];
% p6=bar(log10(COMM));
% hold on
% e=errorbar(log10(COMM), log10(stdCOMM)); %use log10 otherwise communicability is much too large to visualise
% e.LineStyle='none';
% % xlim([0.05 0.6])
% % ylim([2 16])
% xticklabels({[num2str(sizeNetwork) 'node Random Nw'],[num2str(sizeNetwork) 'node Ordered Nw'],'100nw','500nw','1000nw','2000nw'});
% ylabel('Log10 Communicability');
fprintf('Figure 6 Complete \n');
%%
%Betweenness Centrality

f7=figure;
BCRandom=[];
BCOrdered=[];
BCRandomstd=[];
BCOrderedstd=[];
for i = 1:plotNet
    %Circuit Rank:
    BCRandom=[BCRandom [randomOrdered100.AvgBC(i,:)]];
    %     BCOrdered=[BCOrdered ordered100.AvgBC];
    BCRandomstd=[BCRandomstd [randomOrdered100.StdBC(i,:)]];
    %     BCOrderedstd=[BCOrderedstd ordered100.StdBC];
    randomLabel{i}=[num2str(sizeNetwork) ' Random Nw'];
    %     orderedLabel{i}=[num2str(sizeNetwork) ' Ordered Nw'];
end
%BCOrdered BCOrderedstd
BC=[BCRandom  ANN.avgBC cElegans.avgBC AgNW.AvgBC];
stdBC=[BCRandomstd  ANN.stdBC cElegans.stdBC AgNW.StdBC];
p7=bar(BC);
hold on
e=errorbar(BC, stdBC);
e.LineStyle='none';
% xlim([0.05 0.6])
% ylim([2 16])
%orderedLabel{:}
xticklabels([randomLabel{:},'500node Artificial Neural Nw','C. Elegans Nw','100nw','500nw','1000nw','2000nw']);
set(gca, 'XTickLabelRotation', 45)
ylabel('Betweenness Centrality');
hold on


fprintf('Figure 7 Complete \n');



%% SAVE GRAPHS
set(f,'PaperPositionMode','auto');
set(f,'PaperOrientation','landscape');
set(f,'Position',[0 0 1920 1080]);
print(f,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Clustering Coefficient vs Path Length all networks.pdf']);

set(fDeg,'PaperPositionMode','auto');
set(fDeg,'PaperOrientation','landscape');
set(fDeg,'Position',[0 0 1920 1080]);
print(fDeg,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Degree Distribution all networks.pdf']);


set(f2,'PaperPositionMode','auto');
set(f2,'PaperOrientation','landscape');
set(f2,'Position',[0 0 1920 1080]);
print(f2,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Small World Prop all networks.pdf']);

set(f3,'PaperPositionMode','auto');
set(f3,'PaperOrientation','landscape');
set(f3,'Position',[0 0 1920 1080]);
print(f3,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Circuit Rank all networks.pdf']);

set(f4,'PaperPositionMode','auto');
set(f4,'PaperOrientation','landscape');
set(f4,'Position',[0 0 1920 1080]);
print(f4,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Participant Coefficient vs Module z-Score all networks.pdf']);

set(f5,'PaperPositionMode','auto');
set(f5,'PaperOrientation','landscape');
set(f5,'Position',[0 0 1920 1080]);
print(f5,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Guimera PC vs Module z-Score all networks.pdf']);
% print(f6,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Communicability all networks.pdf']);

set(f7,'PaperPositionMode','auto');
set(f7,'PaperOrientation','landscape');
set(f7,'Position',[0 0 1920 1080]);
print(f7,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Betweenness Centrality all networks.pdf']);

end
function plotMulti(sizeNetwork,randomOrdered,randomOrdered100, human, e100Multi, e500Multi, e1000Multi, e2000Multi, AgNW,cElegans,ANN,fig_dir,loop)
% if we are looping through all networks:
f=figure;

while 1
    plotNet=input('Which Network size would you like to plot? 100, 500, 1000, or 2000? \n');
    switch plotNet
        case 100
            plotNet=1;
        case 500
            plotNet=2;
        case 1000
            plotNet=3;
        case 2000
            plotNet=4;
    end
    AgNWidx=length(randomOrdered.GlobalClust(plotNet,:))+3;
    %                 ANN.GlobalClust_Trained ANN.AvgPath_Trained
    y=[randomOrdered.GlobalClust{plotNet,:} ANN.GlobalClust  cElegans.GlobalClust  AgNW.AvgGlobalClust];
    x=[randomOrdered.AvgPath{plotNet,:} ANN.AvgPath  cElegans.AvgPath AgNW.AvgAvgPath];
    p=gscatter(x(1:length([randomOrdered.GlobalClust{plotNet,:}])),y(1:length([randomOrdered.GlobalClust{plotNet,:}])),[],'blue');
    %     set(p, 'linestyle', '-');
    hold on
    p1=gscatter(x(length([randomOrdered.GlobalClust{plotNet,:}])+1:end),y(length([randomOrdered.GlobalClust{plotNet,:}])+1:end),[],'red');
    %     e1=errorbar(x(1), y(1),randomOrdered100.StdPath(plotNet));
    %     e2=errorbar(x(1), y(1),randomOrdered100.StdGlobalClust(plotNet),'horizontal');
    
    % ylim([2 16])
    for m = 1:length(randomOrdered.GlobalClust(plotNet,:)) %for each watts strogatz point
        if m == 1 | m == length(randomOrdered.GlobalClust(plotNet,:))
            label{m}=[num2str(AgNW.NumWires(plotNet)) ' Node WS: B=' num2str((m-1)*0.05)];
        end
    end
    
    text(x,y,{label{:},'500node ANN Untrained','C. Elegans Nw','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
    ylabel('Global Clustering Coefficient');
    xlabel('Global Mean Path Length');
    %     p.MarkerEdgeColor='b';
    %     p(:,1).MarkerEdgeColor='r';
    
    plot2=lower(input('Would you like to plot another network size? \n','s'));
    if plot2=='n'
        break;
    end
end
e3=errorbar(x(AgNWidx:length(x)),y(AgNWidx:length(y)),AgNW.StdGlobalClust,'horizontal');
e4=errorbar(x(AgNWidx:length(x)),y(AgNWidx:length(y)),AgNW.StdAvgPath);

%     e1=errorbar(x(7:end), y(7:end),AgNW.StdAvgPath);
%     e2=errorbar(x(7:end), y(7:end),AgNW.StdGlobalClust,'horizontal');
e3.LineStyle='none';
e4.LineStyle='none';
e3.LineWidth=1.5;
e4.LineWidth=1.5;
e3.Color='black';
e4.Color='black';
% xlim([0.05 0.6])
fprintf('Figure 1 Complete \n');

%%
%Small World Prop


smallWorldPropRandomOrdered=[];
% smallWorldPropOrdered=[];

for i = 1:plotNet
    %Circuit Rank:
    smallWorldPropRandomOrdered=[smallWorldPropRandomOrdered randomOrdered100.AvgSmallWorldProp(i,:)];
    %     smallWorldPropOrdered=[smallWorldPropOrdered ordered100.AvgSmallWorldProp];
    randomOrderedLabel{i}=[num2str(sizeNetwork(i)) ' Watts-Strogatz Nw'];
    %     orderedLabel{i}=[num2str(sizeNetwork) ' Ordered Nw'];
end
f2=figure;

x=[smallWorldPropRandomOrdered  ANN.SmallWorldProp cElegans.SmallWorldProp AgNW.AvgSmallWorldProp];
p2=bar(x);
hold on
for i = 1:plotNet
    e=errorbar(i,x(i), randomOrdered100.StdSmallWorldProp(i,:));
    %     e2=errorbar(x(2),ordered100.StdSmallWorldProp);
    
    % xlim([0.05 0.6])
    % ylim([2 16])
    e.Color='black';
    e.LineWidth=1.5;
    hold on
    
end
e3=errorbar((plotNet+3):length(x),x(plotNet+3:end),AgNW.StdSmallWorldProp);
e3.LineStyle='none';
e3.Color='black';
e3.LineWidth=2;
xticklabels([randomOrderedLabel(1:plotNet),'500node Artificial Neural Nw','C. Elegans Nw','100nw','500nw','1000nw','2000nw']);

set(gca, 'XTickLabelRotation', 45)

ylabel('Small World Prop');


fprintf('Figure 2 Complete \n');

%%
f3=figure;
circuitRankRandomOrdered=[];
% circuitRankOrdered=[];
for i = 1:plotNet
    %Circuit Rank:
    circuitRankRandomOrdered=[circuitRankRandomOrdered randomOrdered.CircuitRank{i}];
    %     circuitRankOrdered=[circuitRankOrdered ordered(1).CircuitRank];
end

circuitRank=[circuitRankRandomOrdered ANN.CircuitRank AgNW.AvgCircuitRank];
p3=bar(circuitRank);
hold on
e5=errorbar((plotNet+2):length(circuitRank),circuitRank(plotNet+2:end),AgNW.StdCircuitRank);
e5.LineStyle='none';
e5.LineWidth=1.5;
e5.Color='black';
xticklabels([randomOrderedLabel(1:plotNet) {'500node Artificial Neural Nw','100nw', '500nw', '1000nw','2000nw'}]);
set(gca, 'XTickLabelRotation', 45)

ylabel('Circuit Rank');
hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles
for i=1:length(p3)  % iterate over number of bar objects
    hT=[hT text(p3(i).XData+p3(i).XOffset,p3(i).YData,num2str(p3(i).YData.'), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end
hold on
%Average Participation Coefficient & Module z-Score
fprintf('Figure 3 Complete \n');
%%
f4=figure;

PRandom=[];
% POrdered=[];
MZRandom=[];
% MZOrdered=[];
stdPRandom=[];
% stdPOrdered=[];
stdMZRandom=[];
% stdMZOrdered=[];

for i = 1:plotNet
    %Circuit Rank:
    PRandom=[PRandom mean(randomOrdered100.AvgPCoeff(i,:))];
    %     POrdered=[POrdered ordered100.AvgAvgPCoeff];
    stdPRandom=[stdPRandom std(randomOrdered100.AvgPCoeff(i,:))];
    %     stdPOrdered=[stdPOrdered ordered100.StdAvgPCoeff];
    MZRandom=[MZRandom mean(randomOrdered100.AvgMZ(i,:))];
    %     MZOrdered=[MZOrdered ordered100.AvgAvgMZ];
    stdMZRandom=[stdMZRandom std(randomOrdered100.AvgMZ(i,:))];
    %     stdMZOrdered=[stdMZOrdered ordered100.StdAvgMZ];
end

PCoeff=[PRandom ANN.avgP cElegans.avgP human.PLocalHubs human.PLocalHubs human.PConnectorHubs human.PConnectorHubs AgNW.AvgPCoeff];
stdPCoeff=[PRandom ANN.avgP cElegans.avgP, 0, 0, 0, 0, 0, AgNW.StdPCoeff];
stdMZ=[stdMZRandom ANN.stdMZ cElegans.stdMZ,0,0,0,0,0, AgNW.AvgMZ];
MZ=[MZRandom ANN.avgMZ cElegans.avgMZ human.MZHubs human.MZNonHubs human.MZHubs human.MZNonHubs  AgNW.StdMZ];

p4=gscatter(PCoeff,MZ);
% e=errorbar(MZ, stdMZ);
% e2=errorbar(PCoeff,stdPCoeff);

%High PCoeff = Hubs / Central areas (Power et al., 2013)
text(PCoeff,MZ,[randomOrderedLabel(1:plotNet), '500node Artificial Neural Nw', 'C. Elegans Nw', 'Human Connector Local Provincial Hub','Human Local Peripheral Node','Human Connector Hub','Human Satellite Connector', '100nw Avg', '500nw Avg', '1000nw Avgs','2000nw Avg'],'NorthWest');
xlabel('Average Participant Coefficient Coefficient');
ylabel('Average Module z-Score');
p4.MarkerEdgeColor='b';
p4(:,1).MarkerEdgeColor='r';
p4.LineWidth=1.5;
hold on

fprintf('Figure 4 Complete \n');

%% Plot Guimera & Amaral rectangles:
f5=figure;
for i=1:plotNet
    subplot(2,2,i)
    if i == 1
        network{i}=e100Multi;
    elseif i==2
        network{i}=e500Multi;
    elseif i==3
        network{i}=e1000Multi;
    elseif i==4
        network{i}=e2000Multi;
    end
    guimera(network{i},randomOrdered100,i,sizeNetwork(i),AgNW); %change network here
    fprintf(['Figure 5 part ' num2str(i) ' Complete \n']);
end

fprintf('Figure 5 Complete \n');
%%
% %%
% %Betweenness Centrality
%
% f7=figure;
% BCRandom=[];
% % BCOrdered=[];
% BCRandomstd=[];
% % BCOrderedstd=[];
% for i = 1:plotNet
%     %Circuit Rank:
%     BCRandom=[BCRandom randomOrdered100.AvgBC(i,:)];
%     %     BCOrdered=[BCOrdered ordered100.AvgBC];
%     BCRandomstd=[BCRandomstd randomOrdered100.StdBC(i,:)];
%     %     BCOrderedstd=[BCOrderedstd ordered100.StdBC];
%     randomOrderedLabel{i}=[num2str(sizeNetwork(i)) ' Watts-Strogatz Nw'];
%     %     orderedLabel{i}=[num2str(sizeNetwork(i)) ' Ordered Nw'];
% end
% %
% BC=[BCRandom ANN.avgBC cElegans.avgBC AgNW.AvgBC];
% stdBC=[BCRandomstd ANN.stdBC cElegans.stdBC AgNW.StdBC];
% p7=bar(BC);
% hold on
% e=errorbar(BC, stdBC);
% e.LineStyle='none';
% % xlim([0.05 0.6])
% % ylim([2 16])
% %
%
% xticklabels({randomOrderedLabel{:},'500node Artificial Neural Nw','C. Elegans Nw','100nw','500nw','1000nw','2000nw'});
% set(gca, 'XTickLabelRotation', 45)
% ylabel('Betweenness Centrality');
% hold on
%

% fprintf('Figure 7 Complete \n');
% set(f,'PaperPositionMode','auto');
% set(f,'PaperOrientation','landscape');
% set(f,'Position',[0 0 1920 1080]);
% print(f,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Clust vs PL all networks w Variance.pdf']);
% set(f2,'PaperPositionMode','auto');
% set(f2,'PaperOrientation','landscape');
% set(f2,'Position',[0 0 1920 1080]);
% print(f2,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Small World Prop all networks w Variance.pdf']);
% set(f3,'PaperPositionMode','auto');
% set(f3,'PaperOrientation','landscape');
% set(f3,'Position',[0 0 1920 1080]);
% print(f3,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Circuit Rank all networks w Variance.pdf']);
% set(f4,'PaperPositionMode','auto');
% set(f4,'PaperOrientation','landscape');
% set(f4,'Position',[0 0 1920 1080]);
% print(f4,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'PC vs MZ all networks w Variance.pdf']);
set(f5,'PaperPositionMode','auto');
set(f5,'PaperOrientation','landscape');
set(f5,'Position',[0 0 1920 1080]);
print(f5,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'PC vs MZ Guimera-Amaral all networks.pdf']);
% print(f7,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Betweenness Centrality all networks with Variance.pdf']);
end

function guimera(network, randomOrdered100,plotNet,sizeNetwork,varargin)

dbstop if error
% Set up rectangles:
xR1=[0 0.025 0.025 0];
xR2=[0.025 0.625 0.625 0.025];
xR3=[0.625 0.8 0.8 0.625];
xR4=[0.8 1 1 0.8];
xR5=[0 0.36 0.36 0];
xR6=[0.36 0.75 0.75 0.36];
xR7=[0.75 1 1 0.75];
yR1234=[-3 -3 2.5 2.5];
yR567=[2.5 2.5 8 8];

R1=patch(xR1,yR1234,'black','LineStyle','none','FaceAlpha',0.2);
hold on
R2=patch(xR2,yR1234,'red','LineStyle','none','FaceAlpha',0.2);
R3=patch(xR3,yR1234,'green','LineStyle','none','FaceAlpha',0.2);
R4=patch(xR4,yR1234,'blue','LineStyle','none','FaceAlpha',0.2);
R5=patch(xR5,yR567,'yellow','LineStyle','none','FaceAlpha',0.2);
R6=patch(xR6,yR567,[255 204 153]./255,'LineStyle','none','FaceAlpha',0.2);
R7=patch(xR7,yR567,[0.7 0.7 0.7],'LineStyle','none','FaceAlpha',0.2);

PCoeff=zeros(length(network.Explore),sizeNetwork);
MZ=zeros(length(network.Explore),sizeNetwork);

% Plot Guimera & Amaral Participant Coefficient & Module z-Score:
if ~isempty(varargin)
    for j = 1:length(network.Explore)
        A=PCoeff(j,:);
        B=[network.Explore{j}.GraphTheory.P]';
        A(padarray(true(size(B)),size(A)-size(B),'post')) = B;
        A2=MZ(j,:);
        B2=[network.Explore{j}.GraphTheory.MZ]';
        A2(padarray(true(size(B2)),size(A2)-size(B2),'post')) = B2;
        
        PCoeff(j,:)=A;
        MZ(j,:)=A2;
        
    end
    
    avgP=mean(PCoeff);
    stdP=std(PCoeff);
    avgMZ=mean(MZ);
    stdMZ=std(MZ);
%     avgMZ(avgMZ==0)=[];
%     avgP(avgP==0)=[];
    PCoeffRandom=randomOrdered100.AvgPCoeff(plotNet,:);
    % PCoeffOrdered= ordered100.AvgPCoeff;
    MZRandom=randomOrdered100.AvgMZ(plotNet,:);
    % MZOrdered=ordered100.AvgMZ;
    
else
    PCoeff=[network.Explore.GraphTheory.P];
    MZ=[network.Explore.GraphTheory.MZ];
    PCoeffRandom=randomOrdered100.AvgPCoeff;
    % PCoeffOrdered= ordered100.AvgPCoeff;
    MZRandom=randomOrdered100.AvgMZ;
    % MZOrdered=ordered100.AvgMZ;
end
if ~isempty(varargin)
    for i = 1:length(avgP)
        if avgMZ(i)<2.5 && avgP(i)>=0 && avgP(i)<0.025
            RegionsMap(i)=1;
        elseif avgMZ(i)<2.5 && avgP(i)>=0.025 && avgP(i)<=0.625
            RegionsMap(i)=2;
        elseif avgMZ(i)<2. && avgP(i) >0.625 && avgP(i) <=0.8
            RegionsMap(i)=3;
        elseif avgMZ(i)<2.5 && avgP(i) >0.8
            RegionsMap(i)=4;
        elseif avgMZ(i)>=2.5 && avgP(i)>=0 && avgP(i)<=0.36
            RegionsMap(i)=5;
        elseif avgMZ(i)>=2.5 && avgP(i)>0.36 && avgP(i)<=0.75
            RegionsMap(i)=6;
        elseif avgMZ(i)>=2.5 && avgP(i) >0.75
            RegionsMap(i)=7;
        end
    end
    
    hold on;
    
       %Find the max std value out of all the values:
    [c,idxMax]=max(avgP+stdP);
    [cmin,idxMin]=min(avgP-stdP);
    [c2,idx2Max]=max(avgMZ+stdMZ);
    [c2min,idx2Min]=min(avgMZ-stdMZ);
    
    l=line([cmin, c],[(c2min+c2)/2,(c2min+c2)/2]);
     l2=line([(cmin+c)/2,(cmin+c)/2],[c2min, c2]);
     l.Color='black';
     l2.Color='black';
     l.LineWidth=1.5;
     l2.LineWidth=1.5;
        hold on;

    for i = 1:length(avgP)
        switch RegionsMap(i)
            case 1
                r1=scatter(avgP(i),avgMZ(i),'black','filled');
                r1.MarkerEdgeColor='black';
            case 2
                r2=scatter(avgP(i),avgMZ(i),'red','filled');
                r2.MarkerEdgeColor='black';
            case 3
                r3=scatter(avgP(i),avgMZ(i),'green','filled');
                r3.MarkerEdgeColor='black';
            case 4
                r4=scatter(avgP(i),avgMZ(i),'blue','filled');
                r4.MarkerEdgeColor='black';
            case 5
                r5=scatter(avgP(i),avgMZ(i),'yellow','filled');
                r5.MarkerEdgeColor='black';
            case 6
                r6=scatter(avgP(i),avgMZ(i),[],[255 204 153]/255,'filled');
                r6.MarkerEdgeColor='black';
            case 7
                r7=scatter(avgP(i),avgMZ(i),[],[0.7 0.7 0.7],'filled');
                r7.MarkerEdgeColor='black';
        end
    end
 
%                     e1h=errorbar(idxMax,[cmin:c],'horizontal');
%                 e1=errorbar(avgP,avgMZ,stdMZ);
else
    for i = 1:length(PCoeff)
        if MZ(i)<2.5 && PCoeff(i)>=0 && PCoeff(i)<0.025
            RegionsMap(i)=1;
        elseif MZ(i)<2.5 && PCoeff(i)>=0.025 && PCoeff(i)<=0.625
            RegionsMap(i)=2;
        elseif MZ(i)<2. && PCoeff(i) >0.625 && PCoeff(i) <=0.8
            RegionsMap(i)=3;
        elseif MZ(i)<2.5 && PCoeff(i) >0.8
            RegionsMap(i)=4;
        elseif MZ(i)>=2.5 && PCoeff(i)>=0 && PCoeff(i)<=0.36
            RegionsMap(i)=5;
        elseif MZ(i)>=2.5 && PCoeff(i)>0.36 && PCoeff(i)<=0.75
            RegionsMap(i)=6;
        elseif MZ(i)>=2.5 && PCoeff(i) >0.75
            RegionsMap(i)=7;
        end
    end
    hold on;
    for i = 1:length(PCoeff)
        switch RegionsMap(i)
            case 1
                r1=scatter(PCoeff(i),MZ(i),'black','filled');
                r1.MarkerEdgeColor='black';
            case 2
                r2=scatter(PCoeff(i),MZ(i),'red','filled');
                r2.MarkerEdgeColor='black';
            case 3
                r3=scatter(PCoeff(i),MZ(i),'green','filled');
                r3.MarkerEdgeColor='black';
            case 4
                r4=scatter(PCoeff(i),MZ(i),'blue','filled');
                r4.MarkerEdgeColor='black';
            case 5
                r5=scatter(PCoeff(i),MZ(i),'yellow','filled');
                r5.MarkerEdgeColor='black';
            case 6
                r6=scatter(PCoeff(i),MZ(i),[],[255 204 153]/255,'filled');
                r6.MarkerEdgeColor='black';
            case 7
                r7=scatter(PCoeff(i),MZ(i),[],[0.7 0.7 0.7],'filled');
                r7.MarkerEdgeColor='black';
        end
    end
end

ylim([-3 8]);
xlim([0,1]);
xlabel('Participant Coefficient')
ylabel('Within-Module Degree z-Score');
title(['Avg Across 11 ' num2str(sizeNetwork) 'nw Networks']);
hold on
end