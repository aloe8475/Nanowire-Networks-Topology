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
% Betweenness centrality measures.
%-------------------------------------------------------------------------

dbstop if error
close all


loadData=lower(input('Would you like to load cross-network data? \n','s'));
if loadData=='y'
    loadDataMulti=lower(input('Would you like to compare multiple iterations of a network? \n','s'));
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
                    explore_location_origin='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\';
                    explore_location=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\MultiNetwork Analysis\' num2str(currMultiNet) 'nw Alternate NWs\'];
                    savePath=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Multi-Network Data\'];
                    fig_dir=['C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\' num2str(currMultiNet) 'nw Alternate NWs\'];
                end
            case '' %if on linux
                explore_location='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/';
                savePath='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/Multi-Network Data/';
                fig_dir='/suphys/aloe8475/Documents/CODE/Data/Figures/Explore Analysis/Cross-Network Explore/Graph Theory/';
                
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
            e100=load([explore_location 'Zdenka_Net_Sx20_NoW109_0930-2019_110329_Length_6_Disp_0_Sim_1_Source_99_Drain_36_Explore_Timestamp_200.mat']);
            e500=load([explore_location 'Zdenka_Net_Sx20_NoW500_0330-2019_111659_Length_6_Disp_0_Sim_1_Source_498_Drain_435_Explore_Timestamp_200.mat']);
            e1000=load([explore_location 'Zdenka_Net_Sx20_NoW1000_0606-2019_113353_Length_6_Disp_0_Sim_1_Source_997_Drain_408_Explore_Timestamp_200.mat']);
            e2000=load([explore_location 'Zdenka_Net_Sx20_NoW2000_0618-2019_125103_Length_6_Disp_0_Sim_1_Source_1999_Drain_935_Explore_Timestamp_200.mat']);
            
            fprintf('Data Loaded \n');
            
            cd(explore_location);
            fprintf(['Loading Multi Network Data from Network size ' num2str(currMultiNet) '\n']);
            dinfo = dir('Adrian*.mat');
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
        [Net, random, ordered,network]=randomOrdered(savePath,currentLocation,AgNW);%e100,e500,e1000,e2000);
    end
  
    %Plot graphs
    plotAll(Net,random,ordered, human, e100, e500, e1000, e2000, AgNW, network,cElegans,ANN,fig_dir)
    
else % if we are looping through
    loop = 1;
    ANN=AI(e500,savePath,loop);
    cElegans=cElegansFun();
    human=humanFun();
    if loadData == 'y'
        [Net, random, ordered,network]=randomOrdered(savePath,currentLocation,e100,e500,e1000,e2000);
    end
    AgNW=AgNWFunMulti(e100Multi, e500Multi, e1000Multi, e2000Multi);    %function for looping
    
    plotMulti(Net,random,ordered, human, e100Multi, e500Multi, e1000Multi, e2000Multi, AgNW, network,cElegans,ANN,fig_dir)
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
function [Net, random, ordered,network]=randomOrdered(savePath,currentLocation,AgNW)
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
%             Net(i).sizeNetwork=100;
%         case 500
%             i=2;
%             network{i}=e500;
%             
%             Net(i).sizeNetwork=500;
%         case 1000
%             i=3;
%             network{i}=e1000;
%             
%             Net(i).sizeNetwork=1000;
%         case 2000
%             i=4;
%             network{i}=e2000;
%             
%             Net(i).sizeNetwork=2000;
%     end
for k = 1:length(AgNW.NumWires)
    Net(k).sizeNetwork=AgNW.NumWires(k);
end 
    if createNewRand=='n'
        if exist([loadPath 'Random_Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'], 'file') == 2 | exist([loadPath 'Random_Ordered_Graphs_ALL_networks.mat'], 'file') == 2 %&& exist([loadPath 'Random_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'],'file') == 2 %2 because .mat file
            while 1
                [a b]=find([Net.sizeNetwork]==Net(i).sizeNetwork);
                if sum(a)>1 %if we load the same network twice
                    waitfor(msgbox('Error: Cannot load same network twice. Please select a different Network'));
                else
                    fprintf(['Loading Ordered and Random Graphs (' num2str(Net(i).sizeNetwork) 'nodes)... \n \n']);
%                     load([loadPath 'Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
                    load([loadPath 'Random_Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
                    break;
                end
            end
        else
            fprintf('Ordered and Random Graphs have not been created yet \n');
            fprintf('Creating New Graphs \n');
%             [random, ordered, randomOrdered100, ordered100, Parameters]=createRandom_Ordered_Graphs(Net(i).sizeNetwork,AgNW);
            [random, ordered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW);

        end
    else
        if exist([loadPath 'Random_Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'], 'file') == 2 | exist([loadPath 'Random_Ordered_Graphs_ALL_networks.mat'], 'file') == 2 %%&& exist([loadPath 'Random_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'],'file') == 2 %2 because .mat file
            overwrite=lower(input(['Graphs (' num2str(Net(i).sizeNetwork) 'nodes) have already been created, would you like to overwrite? \n'],'s'));
            if overwrite =='n'
                fprintf(['Loading Ordered and Random Graphs (' num2str(Net(i).sizeNetwork) 'nodes) \n']);
%                 load([loadPath 'Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
                load([loadPath 'Random_Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
            else
                fprintf('Creating New Graphs \n');
%                 [random, ordered, randomOrdered100, ordered100, Parameters]=createRandom_Ordered_Graphs(Net(i).sizeNetwork,AgNW);
                [random, ordered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW);

                fprintf('New Graphs Created \n');
            end
        else
            fprintf('Creating New Graphs \n');
            [random, ordered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW);
            fprintf('New Graphs Created \n');
        end
    end
    %save in struct so we can load multiple networks at once:
    if Net(i).sizeNetwork == 100 && i == 1
        Net(i).randomOrdered100=randomOrdered100;
        Net(i).randomOrdered=randomOrdered;
%         Net(i).ordered100=ordered100;
%         Net(i).ordered=ordered;
    elseif Net(i).sizeNetwork == 500 && i == 2
        Net(i).randomOrdered100=randomOrdered100;
        Net(i).randomOrdered=randomOrdered;
%         Net(i).ordered100=ordered100;
%         Net(i).ordered=ordered;
    elseif Net(i).sizeNetwork == 1000 && i == 3
        Net(i).randomOrdered100=randomOrdered100;
        Net(i).randomOrdered=randomOrdered;
%         Net(i).ordered100=ordered100;
%         Net(i).ordered=ordered;
    elseif Net(i).sizeNetwork == 2000 && i == 4
        Net(i).randomOrdered100=randomOrdered100;
        Net(i).randomOrdered=randomOrdered;
%         Net(i).ordered100=ordered100;
%         Net(i).ordered=ordered;
    else
        fprintf('Error in Network sizes');
    end
    load2=lower(input('Would you like to load another network? \n','s'));
    if load2=='n'
        break;
    end
end

fprintf('Graphs Loaded... \n Creating Plots... \n');

%% Random & Ordered
%Number of Nodes
for z=1:length(Net)
    while isempty(Net(z).sizeNetwork) %if we only load 500, 1000 or 2000, we need to increase our counter to skip over the empty arrays.
        z=z+1;
    end
    for i =1:100
        Net(z).randomOrdered(i).numNodes=height(Net(z).randomOrdered(i).Graph.Nodes);
%         Net(z).ordered(i).numNodes=height(Net(z).ordered(i).Graph.Nodes);
        Net(z).randomOrdered(i).numEdges=height(Net(z).randomOrdered(i).Graph.Edges);
%         Net(z).ordered(i).numEdges=height(Net(z).ordered(i).Graph.Edges);
    end
    Net(z).randomOrdered100.AvgNodes=mean([Net(z).randomOrdered(:).numNodes]);
%     Net(z).ordered100.AvgNodes=mean([Net(z).ordered(:).numNodes]);
    
    %Number of Edges
    Net(z).randomOrdered100.AvgEdges=mean([Net(z).randomOrdered(:).numEdges]);
%     Net(z).ordered100.AvgEdges=mean([Net(z).ordered(:).numEdges]);
    
    % participation coefficient (mean)
    
    % avg and std pcoeff across 100 bootraps per node:
    Net(z).randomOrdered100.AvgPCoeff=mean([Net(z).randomOrdered(:).P],2);
%     Net(z).ordered100.AvgPCoeff=mean([Net(z).ordered(:).P],2);
    Net(z).randomOrdered100.StdPCoeff=std([Net(z).randomOrdered(:).P],[],2);
%     Net(z).ordered100.StdPCoeff=std([Net(z).ordered(:).P],[],2);

    % Avg & std PCoeff (across 100 bootstraps) and across nodes:
    Net(z).randomOrdered100.AvgAvgPCoeff=mean(mean([Net(z).randomOrdered(:).P]),2);
%     Net(z).randomOrdered100.StdAvgPCoeff=std(mean([Net(z).random(:).P]),[],2);
    Net(z).ordered100.AvgAvgPCoeff=mean(mean([Net(z).randomOrdered(:).P]),2);
%     Net(z).ordered100.StdAvgPCoeff=std(mean([Net(z).ordered(:).P]),[],2);
    %Module z-Score
    % avg and std mz across 100 bootraps per node:
    Net(z).randomOrdered100.AvgMZ=mean([Net(z).randomOrdered(:).MZ],2);
%     Net(z).ordered100.AvgMZ=mean([Net(z).ordered(:).MZ],2);
    Net(z).randomOrdered100.StdMZ=std([Net(z).randomOrdered(:).MZ],[],2);
%     Net(z).ordered100.StdMZ=std([Net(z).ordered(:).MZ],[],2);
    
    %avg and std mz across 100 bootstraps and across nodes:
    Net(z).randomOrdered100.AvgAvgMZ=mean(mean([Net(z).randomOrdered(:).MZ]),2);
    Net(z).randomOrdered100.StdAvgMZ=std(mean([Net(z).randomOrdered(:).MZ]),[],2);
%     Net(z).ordered100.AvgAvgMZ=mean(mean([Net(z).ordered(:).MZ]),2);
%     Net(z).ordered100.StdAvgMZ=std(mean([Net(z).ordered(:).MZ]),[],2);
    %Small World Prop
    Net(z).randomOrdered100.AvgSmallWorldProp=mean([Net(z).randomOrdered(:).SmallWorldProp],2);
    Net(z).randomOrdered100.StdSmallWorldProp=std([Net(z).randomOrdered(:).SmallWorldProp],[],2);
%     Net(z).ordered100.AvgSmallWorldProp=mean([Net(z).ordered(:).SmallWorldProp],2);
%     Net(z).ordered100.StdSmallWorldProp=std([Net(z).ordered(:).SmallWorldProp],[],2);
    
    % communicability
    
    % complexity (from the Sporns, Tononi and Edelman paper I sent through the other day -- still waiting on code from Olaf, but will send through when I get it).
    
    % betweeness_centrality
    
    %Average Normalised Betweenness (BC/[(N-1)(N-2)])
    %     randomOrdered100.AvgNormBC=mean([random(:).normBC]);
    %     randomOrdered100.StdNormBC=std([random(:).normBC]);
    %     ordered100.AvgNormBC=mean([ordered(:).normBC]);
    %     ordered100.StdNormBC=std([ordered(:).normBC]);
end
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
e100.AvgDEG=mean([Degree]);
e100.StdDEG=std([Degree])./sqrt(length([Degree]));
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
e500.AvgDEG=mean([Degree]);
e500.StdDEG=std([Degree])./sqrt(length([Degree]));
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
e1000.AvgDEG=mean([Degree]);
e1000.StdDEG=std([Degree])./sqrt(length([Degree]));
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
e2000.AvgDEG=mean([Degree]);
e2000.StdDEG=std([Degree])./sqrt(length([Degree]));
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
AgNW.AvgDegree=[e100.AvgDegree e500.AvgDegree e1000.AvgDegree e2000.AvgDegree];
AgNW.StdDegree=[e100.StdDegree e500.StdDegree e1000.StdDegree e2000.StdDegree];
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

end

%% Plot:
function plotAll(Net, random, ordered, human, e100, e500, e1000, e2000, AgNW,network,cElegans,ANN,fig_dir)

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
        y=[Net(plotNet).randomOrdered100.AvgGlobalClust ANN.GlobalClust ANN.GlobalClust_Trained cElegans.GlobalClust human.GlobalClust  e100.Explore.GraphTheory.GlobalClust e500.Explore.GraphTheory.GlobalClust e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
                %Net(plotNet).ordered100.AvgGlobalClust
                %Net(plotNet).ordered100.AvgPath
        x=[Net(plotNet).randomOrdered100.AvgPath ANN.AvgPath ANN.AvgPath_Trained cElegans.AvgPath human.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
        %                p=gscatter(x,y);
        %     hold on
        for i=1:length(x)
            %           s(i+1)=scatter(nanmean(categories.Clust{i}),categories.PathLength{i,j},[],clrs2{1}(i,:));
            s(i)=scatter(x(:,i),y(:,i),[],clrs2{1}(i,:));
            hold on
            e=errorbar(x(1), y(1),Net(plotNet).randomOrdered100.StdPath);
            e2=errorbar(x(1), y(1),Net(plotNet).randomOrdered100.StdGlobalClust,'horizontal');
%             e3=errorbar(x(5), y(5),Net(plotNet).ordered100.StdPath);
%             e4=errorbar(x(5), y(5),Net(plotNet).ordered100.StdGlobalClust,'horizontal');
            e.Color=clrs2{1}(i,:);
            e2.Color=clrs2{1}(i,:);
%             e3.Color=clrs2{1}(i,:);
%             e4.Color=clrs2{1}(i,:);
        end
                text(x,y,{[num2str(Net(plotNet).sizeNetwork) 'node Watts-Strogats'],'500node ANN Untrained','500node ANN Trained','C. Elegans Nw', 'Human Nw (mm-scale regions)','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')

        %     hold on
        %     e=errorbar(x(1), y(1),Net(plotNet).randomOrdered100.StdPath);
        %     e2=errorbar(x(1), y(1),Net(plotNet).randomOrdered100.StdGlobalClust,'horizontal');

        %[num2str(Net(plotNet).sizeNetwork) 'node Ordered Nw']
        
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

% %Small World Log
% f1=figure;
% logx=log10(x);
% logy=log10(y);
% p1=gscatter(logx,y);
% hold on
% hold on
% e=errorbar(logx(1), y(1),Net(plotNet).randomOrdered100.StdPath);
% e2=errorbar(logx(1), y(1),Net(plotNet).randomOrdered100.StdGlobalClust);
% errorbar(logx(3), y(3),Net(plotNet).ordered100.StdPath);
% errorbar(logx(3), y(3),Net(plotNet).ordered100.StdGlobalClust);
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
    smallWorldPropOrderedRandom=[smallWorldPropOrderedRandom Net(i).randomOrdered100.AvgSmallWorldProp];
%     smallWorldPropOrdered=[smallWorldPropOrdered Net(i).ordered100.AvgSmallWorldProp];
    randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Watts-Strogatz Nw'];
%     orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
end

%Small World Prop
f2=figure;
%smallWorldPropOrdered
x=[smallWorldPropOrderedRandom  ANN.SmallWorldProp cElegans.SmallWorldProp AgNW.SmallWorldProp];
p2=bar(x);
hold on
for i = 1:plotNet
    e=errorbar(x(1), Net(i).randomOrdered100.StdSmallWorldProp);
%     e2=errorbar(x(2),Net(i).ordered100.StdSmallWorldProp);
    % xlim([0.05 0.6])
    % ylim([2 16])
    hold on
end
%orderedLabel
xticklabels([randomLabel,'500node ANN','C. Elegans Nw','100nw','500nw','1000nw','2000nw']);

set(gca, 'XTickLabelRotation', 45)

ylabel('Small World Prop');


fprintf('Figure 2 Complete \n');


f3=figure;
circuitRankRandom=[];
circuitRankOrdered=[];
for i = 1:plotNet
    %Circuit Rank:
    circuitRankRandom=[circuitRankRandom Net(i).randomOrdered(1).CircuitRank];
%     circuitRankOrdered=[circuitRankOrdered Net(i).ordered(1).CircuitRank];
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
    PRandom=[PRandom Net(i).randomOrdered100.AvgAvgPCoeff];
%     POrdered=[POrdered Net(i).ordered100.AvgAvgPCoeff];
    stdPRandom=[stdPRandom Net(i).randomOrdered100.StdAvgPCoeff];
%     stdPOrdered=[stdPOrdered Net(i).ordered100.StdAvgPCoeff];
    MZRandom=[MZRandom Net(i).randomOrdered100.AvgAvgMZ];
%     MZOrdered=[MZOrdered Net(i).ordered100.AvgAvgMZ];
    stdMZRandom=[stdMZRandom Net(i).randomOrdered100.StdAvgMZ];
%     stdMZOrdered=[stdMZOrdered Net(i).ordered100.StdAvgMZ];
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
    guimera(network{i},Net,i); %change network here
    fprintf(['Figure 5 part ' num2str(i) ' Complete \n']);
end

fprintf('Figure 5 Complete \n');

% %Communicability
% f6=figure;
% COMM=[Net(plotNet).randomOrdered100.AvgCOMM(1) Net(plotNet).ordered100.AvgCOMM(1) AgNW.AvgCOMM];
% stdCOMM=[0 0 AgNW.StdCOMM];
% p6=bar(log10(COMM));
% hold on
% e=errorbar(log10(COMM), log10(stdCOMM)); %use log10 otherwise communicability is much too large to visualise
% e.LineStyle='none';
% % xlim([0.05 0.6])
% % ylim([2 16])
% xticklabels({[num2str(Net(plotNet).sizeNetwork) 'node Random Nw'],[num2str(Net(plotNet).sizeNetwork) 'node Ordered Nw'],'100nw','500nw','1000nw','2000nw'});
% ylabel('Log10 Communicability');
fprintf('Figure 6 Complete \n');

%Betweenness Centrality

f7=figure;
BCRandom=[];
BCOrdered=[];
BCRandomstd=[];
BCOrderedstd=[];
for i = 1:plotNet
    %Circuit Rank:
    BCRandom=[BCRandom Net(i).randomOrdered100.AvgBC];
%     BCOrdered=[BCOrdered Net(i).ordered100.AvgBC];
    BCRandomstd=[BCRandomstd Net(i).randomOrdered100.StdBC];
%     BCOrderedstd=[BCOrderedstd Net(i).ordered100.StdBC];
    randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Random Nw'];
%     orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
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
function plotMulti(Net,random,ordered, human, e100Multi, e500Multi, e1000Multi, e2000Multi, AgNW, network,cElegans,ANN,fig_dir,loop)
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
    
    y=[Net(plotNet).random100.AvgGlobalClust ANN.GlobalClust ANN.GlobalClust_Trained cElegans.GlobalClust human.GlobalClust Net(plotNet).ordered100.AvgGlobalClust AgNW.AvgGlobalClust];
    x=[Net(plotNet).random100.AvgPath ANN.AvgPath ANN.AvgPath_Trained cElegans.AvgPath human.AvgPath Net(plotNet).ordered100.AvgPath AgNW.AvgAvgPath];
    p=gscatter(x,y);
    hold on
    errorbar(x(1), y(1),Net(plotNet).random100.StdPath);
    errorbar(x(1), y(1),Net(plotNet).random100.StdGlobalClust,'horizontal');
    errorbar(x(6), y(6),Net(plotNet).ordered100.StdPath);
    errorbar(x(6), y(6),Net(plotNet).ordered100.StdGlobalClust,'horizontal');
    e1=errorbar(x(7:end), y(7:end),AgNW.StdAvgPath);
    e2=errorbar(x(7:end), y(7:end),AgNW.StdGlobalClust,'horizontal');
    e1.LineStyle='none';
    e2.LineStyle='none';
    % xlim([0.05 0.6])
    % ylim([2 16])
    text(x,y,{[num2str(Net(plotNet).sizeNetwork) 'node Random Nw'],'500node ANN Untrained','500node ANN Trained','C. Elegans Nw', 'Human Nw (mm-scale regions)',[num2str(Net(plotNet).sizeNetwork) 'node Ordered Nw'],'100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
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

smallWorldPropRandom=[];
smallWorldPropOrdered=[];

for i = 1:plotNet
    %Circuit Rank:
    smallWorldPropRandom=[smallWorldPropRandom Net(i).random100.AvgSmallWorldProp];
    smallWorldPropOrdered=[smallWorldPropOrdered Net(i).ordered100.AvgSmallWorldProp];
    randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Random Nw'];
    orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
end

%Small World Prop
f2=figure;

x=[smallWorldPropRandom smallWorldPropOrdered  ANN.SmallWorldProp cElegans.SmallWorldProp AgNW.AvgSmallWorldProp];
p2=bar(x);
hold on
for i = 1:plotNet
    e=errorbar(x(1), Net(i).random100.StdSmallWorldProp);
    e2=errorbar(x(2),Net(i).ordered100.StdSmallWorldProp);
    e3=errorbar(x(5:end),AgNW.StdSmallWorldProp);
    % xlim([0.05 0.6])
    % ylim([2 16])
    hold on
end
xticklabels([randomLabel,orderedLabel,'500node Artificial Neural Nw','C. Elegans Nw','100nw','500nw','1000nw','2000nw']);

set(gca, 'XTickLabelRotation', 45)

ylabel('Small World Prop');


fprintf('Figure 2 Complete \n');


f3=figure;
circuitRankRandom=[];
circuitRankOrdered=[];
for i = 1:plotNet
    %Circuit Rank:
    circuitRankRandom=[circuitRankRandom Net(i).random(1).CircuitRank];
    circuitRankOrdered=[circuitRankOrdered Net(i).ordered(1).CircuitRank];
end

circuitRank=[circuitRankRandom circuitRankOrdered ANN.CircuitRank AgNW.AvgCircuitRank];
p3=bar(circuitRank);
errorbar(circuitRank(4:end), AgNW.StdSmallWorldProp);
xticklabels([randomLabel orderedLabel {'500node Artificial Neural Nw','100nw', '500nw', '1000nw','2000nw'}]);
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
    PRandom=[PRandom Net(i).random100.AvgAvgPCoeff];
    POrdered=[POrdered Net(i).ordered100.AvgAvgPCoeff];
    stdPRandom=[stdPRandom Net(i).random100.StdAvgPCoeff];
    stdPOrdered=[stdPOrdered Net(i).ordered100.StdAvgPCoeff];
    MZRandom=[MZRandom Net(i).random100.AvgAvgMZ];
    MZOrdered=[MZOrdered Net(i).ordered100.AvgAvgMZ];
    stdMZRandom=[stdMZRandom Net(i).random100.StdAvgMZ];
    stdMZOrdered=[stdMZOrdered Net(i).ordered100.StdAvgMZ];
end

PCoeff=[PRandom POrdered  ANN.avgP cElegans.avgP human.AvgP human.PLocalHubs human.PLocalHubs human.PConnectorHubs human.PConnectorHubs  mean(AgNW.AvgPCoeff)];
stdPCoeff=[PRandom POrdered  ANN.avgP cElegans.avgP, 0, 0, 0, 0, 0,  std(AgNW.AvgPCoeff)];
stdMZ=[stdMZRandom stdMZOrdered ANN.stdMZ cElegans.stdMZ,0,0,0,0,0, mean(AgNW.AvgMZ)];
MZ=[MZRandom MZOrdered ANN.avgMZ cElegans.avgMZ human.AvgMZ human.MZHubs human.MZNonHubs human.MZHubs human.MZNonHubs  std(AgNW.AvgMZ)];

p4=gscatter(PCoeff,MZ);
% e=errorbar(MZ, stdMZ);
% e2=errorbar(PCoeff,stdPCoeff);

%High PCoeff = Hubs / Central areas (Power et al., 2013)
text(PCoeff,MZ,[randomLabel, orderedLabel, '500node Artificial Neural Nw', 'C. Elegans Nw', 'Human Average', 'Human Connector Local Provincial Hub','Human Local Peripheral Node','Human Connector Hub','Human Satellite Connector', '100nw Avg', '500nw Avg', '1000nw Avgs','2000nw Avg'],'NorthWest');
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
    guimera(network{i},Net,i,AgNW); %change network here
    fprintf(['Figure 5 part ' num2str(i) ' Complete \n']);
end

fprintf('Figure 5 Complete \n');


%Betweenness Centrality

f7=figure;
BCRandom=[];
BCOrdered=[];
BCRandomstd=[];
BCOrderedstd=[];
for i = 1:plotNet
    %Circuit Rank:
    BCRandom=[BCRandom Net(i).random100.AvgBC];
    BCOrdered=[BCOrdered Net(i).ordered100.AvgBC];
    BCRandomstd=[BCRandomstd Net(i).random100.StdBC];
    BCOrderedstd=[BCOrderedstd Net(i).ordered100.StdBC];
    randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Watts-Strogatz Nw'];
    orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
end
% 
BC=[BCRandom BCOrdered ANN.avgBC cElegans.avgBC AgNW.AvgBC];
stdBC=[BCRandomstd BCOrderedstd  ANN.stdBC cElegans.stdBC AgNW.StdBC];
p7=bar(BC);
hold on
e=errorbar(BC, stdBC);
e.LineStyle='none';
% xlim([0.05 0.6])
% ylim([2 16])
%

xticklabels({randomLabel{:},orderedLabel{:},'500node Artificial Neural Nw','C. Elegans Nw','100nw','500nw','1000nw','2000nw'});
set(gca, 'XTickLabelRotation', 45)
ylabel('Betweenness Centrality');
hold on


fprintf('Figure 7 Complete \n');
set(f,'PaperPositionMode','auto');
set(f,'PaperOrientation','landscape');
set(f,'Position',[0 0 1920 1080]);
print(f,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Clustering Coefficient vs Path Length all networks with Variance.pdf']);
set(f2,'PaperPositionMode','auto');
set(f2,'PaperOrientation','landscape');
set(f2,'Position',[0 0 1920 1080]);
print(f2,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Small World Prop all networks with Variance.pdf']);
set(f3,'PaperPositionMode','auto');
set(f3,'PaperOrientation','landscape');
set(f3,'Position',[0 0 1920 1080]);
print(f3,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Circuit Rank all networks with Variance.pdf']);
set(f4,'PaperPositionMode','auto');
set(f4,'PaperOrientation','landscape');
set(f4,'Position',[0 0 1920 1080]);
print(f4,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Participant Coefficient vs Module z-Score all networks with Variance.pdf']);
% print(f6,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Communicability all networks.pdf']);
print(f7,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Betweenness Centrality all networks with Variance.pdf']);
end

function guimera(network, Net,plotNet,varargin)
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

% Plot Guimera & Amaral Participant Coefficient & Module z-Score:

PCoeff=[network.Explore.GraphTheory.P];
MZ=[network.Explore.GraphTheory.MZ];
PCoeffRandom=Net(plotNet).randomOrdered100.AvgPCoeff;
PCoeffOrdered= Net(plotNet).ordered100.AvgPCoeff;
MZRandom=Net(plotNet).randomOrdered100.AvgMZ;
MZOrdered=Net(plotNet).ordered100.AvgMZ;

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
if isempty(varargin)
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
else
    for i = 1:length(varargin)
        switch RegionsMap(i)
            case 1
                r1=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),'black','filled');
                r1.MarkerEdgeColor='black';
            case 2
                r2=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),'red','filled');
                r2.MarkerEdgeColor='black';
            case 3
                r3=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),'green','filled');
                r3.MarkerEdgeColor='black';
            case 4
                r4=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),'blue','filled');
                r4.MarkerEdgeColor='black';
            case 5
                r5=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),'yellow','filled');
                r5.MarkerEdgeColor='black';
            case 6
                r6=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),[],[255 204 153]/255,'filled');
                r6.MarkerEdgeColor='black';
            case 7
                r7=scatter(AgNW.AvgPCoeff(i),AgNW.AvgMZ(i),[],[0.7 0.7 0.7],'filled');
                r7.MarkerEdgeColor='black';
        end
    end
end
ylim([-3 8]);
xlabel('Participant Coefficient')
ylabel('Within-Module Degree z-Score');
title([num2str(Net(plotNet).sizeNetwork) 'nw Network']);
hold on
end