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


close all

currentLocation=pwd;
computer=getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        explore_location='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\';
        savePath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Multi-Network Data\';
        fig_dir='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\';

    case '' %if on linux
        explore_location='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/';
        savePath='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/Multi-Network Data/';
        fig_dir='/suphys/aloe8475/Documents/CODE/Data/Figures/Explore Analysis/Cross-Network Explore/Graph Theory/';

    case 'LAPTOP-S1BV3HR7'
        explore_location='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\';
        savePath='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\Multi-Network Data\';
        fig_dir='D:\alon_\Research\PhD\CODE\Data\Figures\Explore Analysis\Cross-Network Explore\Graph Theory\';
        
        %case '' %--- Add other computer paths (e.g. Mike)
end
cd(explore_location);

loadData=lower(input('Would you like to load the network data? \n','s'));
if loadData=='y'
    %Load three explore analyses:
    fprintf('Loading Data... \n');
    clear e100 e500 e1000 e2000
    e100=load([explore_location 'Adrian_Net_Sx20_NoW100_0325-2019_112338__Sim_1_SourceElectrode_6_DrainElectrode_76_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    e500=load([explore_location 'Adrian_Net_Sx20_NoW500_0330-2019_111659__Sim_1_SourceElectrode_18_DrainElectrode_430_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    e1000=load([explore_location 'Adrian_Net_Sx20_NoW1000_0606-2019_113353__Sim_1_SourceElectrode_32_DrainElectrode_1000_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    e2000=load([explore_location 'Adrian_Net_Sx20_NoW2000_0618-2019_125103__Sim_1_SourceElectrode_158_DrainElectrode_1820_Exploration_Analysis_ Timestamp_400_26-Jun-2019.mat']);
    fprintf('Data Loaded \n');
end


%% Run Functions:
cElegans=cElegansFun();
human=humanFun();
[Net, random, ordered,network]=randomOrdered(savePath,currentLocation,e100,e500,e1000,e2000);
AgNW=AgNWFun(e100, e500, e1000, e2000);

%Plot graphs
plotAll(Net,random,ordered, human, e100, e500, e1000, e2000, AgNW, network,cElegans,fig_dir)


%% TO DO

%500 Node Communicability Graph:
% Net(1).random(1).COMM; %communicability does not change across 100 bootstraps, so we can just use #1
% Net(1).ordered(1).COMM;
% e100.Explore.GraphTheory.COMM;
% e500.Explore.GraphTheory.COMM;
% e1000.Explore.GraphTheory.COMM;
% e2000.Explore.GraphTheory.COMM;


%500 Node Participation Coefficient & Module z-Score
%Plot Graph:

% %   adjacency([e500.Explore.GraphView.Nodes e500.Explore.GraphView.Edges])
%     f4=figure;
%     PCoeff=[random100.AvgPCoeff ordered100.AvgPCoeff, e500.Explore.GraphTheory.P];
%     p3=bar(PCoeff);
%     xticklabels({'500node Random Nw', '500node Ordered Nw', '500nw'});
%     ylabel('Participant Coefficient Coefficient');
%     f5=figure;
%     MZ=[random100.AvgMZ  ordered100.AvgMZ, e500.Explore.GraphTheory.MZ];
%     p4=bar(MZ);
%     xticklabels({'500node Random Nw', '500node Ordered Nw', '500nw'});
%     ylabel('Module z-Score');
%

% Complexity

%% Example AI Graph Analysis (Recurrent Neural Network)
%     function AI()
%     %Create Sample RNN: - NEED TO TALK TO MAC ABOUT THIS
%     %AI.AdjMat=zeros([500 500]);
%     end
%% C-Elegans:
function elegans=cElegansFun()

%Kaiser M, Hilgetag CC (2006) Non-Optimal Component Placement, but Short Processing Paths, due to Long-Distance Projections in Neural Systems. PLoS Computational Biology 7:e95
% Kötter R (2004) Online retrieval, processing, and visualization of primate connectivity data from the CoCoMac database. Neuroinformatics 2: 127-144.
% Choe Y, McCormick BH, Koh W (2004) Network connectivity analysis on the temporally augmented C. elegans web: A pilot study. Society of Neuroscience Abstracts 30:921.9.

cd('../../Data/Organic Networks Connectomes/')
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

%% Human Graph Analysis
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
%% Random and Ordered Graph Analysis
function [Net, random, ordered,network]=randomOrdered(savePath,currentLocation,e100,e500,e1000,e2000)
cd(currentLocation)
i=1;
while 1
    sizeNetwork=input('What Size Random/Ordered Networks would you like to create/load? 100, 500, 1000 or 2000? \n');
    createNewRand=lower(input('Would you like to create new Random and Ordered graphs? (Note this will take 4+ Hours) \n','s'));
    loadPath=savePath;
    % save the data for each network & save the network size in a different
    % variable (sizeNetwork).
    switch sizeNetwork
        case 100
            network=e100;
            i=1;
            Net(i).sizeNetwork=100;
        case 500
            network=e500;
            i=2;
            Net(i).sizeNetwork=500;
        case 1000
            network =e1000;
            i=3;
            Net(i).sizeNetwork=1000;
        case 2000
            network =e2000;
            i=4;
            Net(i).sizeNetwork=2000;
    end
    if createNewRand=='n'
        if exist([loadPath 'Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'], 'file') == 2 && exist([loadPath 'Random_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'],'file') == 2 %2 because .mat file
            fprintf(['Loading Ordered and Random Graphs (' num2str(Net(i).sizeNetwork) 'nodes)... \n \n']);
            load([loadPath 'Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
            load([loadPath 'Random_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
        else
            fprintf('Ordered and Random Graphs have not been created yet \n');
            fprintf('Creating New Graphs \n');
            [random, ordered, random100, ordered100, Parameters]=createRandom_Ordered_Graphs(network,Net(i).sizeNetwork);
        end
    else
        if exist([loadPath 'Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'], 'file') == 2 && exist([loadPath 'Random_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat'],'file') == 2 %2 because .mat file
            overwrite=lower(input(['Graphs (' num2str(Net(i).sizeNetwork) 'nodes) have already been created, would you like to overwrite? \n'],'s'));
            if overwrite =='n'
                fprintf(['Loading Ordered and Random Graphs (' num2str(Net(i).sizeNetwork) 'nodes) \n']);
                load([loadPath 'Ordered_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
                load([loadPath 'Random_Graphs_' num2str(Net(i).sizeNetwork) 'nw.mat']);
            else
                fprintf('Creating New Graphs \n');
                [random, ordered, random100, ordered100, Parameters]=createRandom_Ordered_Graphs(network,Net(i).sizeNetwork);
                fprintf('New Graphs Created \n');
            end
        else
            fprintf('Creating New Graphs \n');
            [random, ordered, random100, ordered100, Parameters]=createRandom_Ordered_Graphs(network,Net(i).sizeNetwork);
            fprintf('New Graphs Created \n');
        end
    end
    %save in struct so we can load multiple networks at once:
    if Net(i).sizeNetwork == 100 && i == 1
        Net(i).random100=random100;
        Net(i).random=random;
        Net(i).ordered100=ordered100;
        Net(i).ordered=ordered;
    elseif Net(i).sizeNetwork == 500 && i == 2
        Net(i).random100=random100;
        Net(i).random=random;
        Net(i).ordered100=ordered100;
        Net(i).ordered=ordered;
    elseif Net(i).sizeNetwork == 1000 && i == 3
        Net(i).random100=random100;
        Net(i).random=random;
        Net(i).ordered100=ordered100;
        Net(i).ordered=ordered;
    elseif Net(i).sizeNetwork == 2000 && i == 4
        Net(i).random100=random100;
        Net(i).random=random;
        Net(i).ordered100=ordered100;
        Net(i).ordered=ordered;
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
        Net(z).random(i).numNodes=height(Net(z).random(i).Graph.Nodes);
        Net(z).ordered(i).numNodes=height(Net(z).ordered(i).Graph.Nodes);
        Net(z).random(i).numEdges=height(Net(z).random(i).Graph.Edges);
        Net(z).ordered(i).numEdges=height(Net(z).ordered(i).Graph.Edges);
    end
    Net(z).random100.AvgNodes=mean([Net(z).random(:).numNodes]);
    Net(z).ordered100.AvgNodes=mean([Net(z).ordered(:).numNodes]);
    
    %Number of Edges
    Net(z).random100.AvgEdges=mean([Net(z).random(:).numEdges]);
    Net(z).ordered100.AvgEdges=mean([Net(z).ordered(:).numEdges]);
    
    % participation coefficient (mean)
    
    % avg and std pcoeff across 100 bootraps per node:
    Net(z).random100.AvgPCoeff=mean([Net(z).random(:).P],2);
    Net(z).ordered100.AvgPCoeff=mean([Net(z).ordered(:).P],2);
    Net(z).random100.StdPCoeff=std([Net(z).random(:).P],[],2);
    Net(z).ordered100.StdPCoeff=std([Net(z).ordered(:).P],[],2);
    % Avg & std PCoeff (across 100 bootstraps) and across nodes:
    Net(z).random100.AvgAvgPCoeff=mean(mean([Net(z).random(:).P]),2);
    Net(z).random100.StdAvgPCoeff=std(mean([Net(z).random(:).P]),[],2);
    Net(z).ordered100.AvgAvgPCoeff=mean(mean([Net(z).ordered(:).P]),2);
    Net(z).ordered100.StdAvgPCoeff=std(mean([Net(z).ordered(:).P]),[],2);
    %Module z-Score
    % avg and std mz across 100 bootraps per node:
    Net(z).random100.AvgMZ=mean([Net(z).random(:).MZ],2);
    Net(z).ordered100.AvgMZ=mean([Net(z).ordered(:).MZ],2);
    Net(z).random100.StdMZ=std([Net(z).random(:).MZ],[],2);
    Net(z).ordered100.StdMZ=std([Net(z).ordered(:).MZ],[],2);
    
    %avg and std mz across 100 bootstraps and across nodes:
    Net(z).random100.AvgAvgMZ=mean(mean([Net(z).random(:).MZ]),2);
    Net(z).random100.StdAvgMZ=std(mean([Net(z).random(:).MZ]),[],2);
    Net(z).ordered100.AvgAvgMZ=mean(mean([Net(z).ordered(:).MZ]),2);
    Net(z).ordered100.StdAvgMZ=std(mean([Net(z).ordered(:).MZ]),[],2);
    %Small World Prop
    Net(z).random100.AvgSmallWorldProp=mean([Net(z).random(:).SmallWorldProp],2);
    Net(z).random100.StdSmallWorldProp=std([Net(z).random(:).SmallWorldProp],[],2);
    Net(z).ordered100.AvgSmallWorldProp=mean([Net(z).ordered(:).SmallWorldProp],2);
    Net(z).ordered100.StdSmallWorldProp=std([Net(z).ordered(:).SmallWorldProp],[],2);
    
    % communicability
    
    % complexity (from the Sporns, Tononi and Edelman paper I sent through the other day -- still waiting on code from Olaf, but will send through when I get it).
    
    % betweeness_centrality
    
    %Average Normalised Betweenness (BC/[(N-1)(N-2)])
    %     random100.AvgNormBC=mean([random(:).normBC]);
    %     random100.StdNormBC=std([random(:).normBC]);
    %     ordered100.AvgNormBC=mean([ordered(:).normBC]);
    %     ordered100.StdNormBC=std([ordered(:).normBC]);
end
end

function AgNW=AgNWFun(e100, e500, e1000, e2000)
%% AgNW
%Circuit Rank
AgNW.CircuitRank=[e100.Explore.GraphTheory.CircuitRank e500.Explore.GraphTheory.CircuitRank e1000.Explore.GraphTheory.CircuitRank e2000.Explore.GraphTheory.CircuitRank];
AgNW.GlobalClust=[e100.Explore.GraphTheory.GlobalClust, e500.Explore.GraphTheory.GlobalClust, e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
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
end
%% Plot:
function plotAll(Net, random, ordered, human, e100, e500, e1000, e2000, AgNW,network,cElegans,fig_dir)
% Small World Analysis
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
    y=[Net(plotNet).random100.AvgGlobalClust cElegans.GlobalClust human.GlobalClust Net(plotNet).ordered100.AvgGlobalClust e100.Explore.GraphTheory.GlobalClust e500.Explore.GraphTheory.GlobalClust e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
    x=[Net(plotNet).random100.AvgPath cElegans.AvgPath human.AvgPath Net(plotNet).ordered100.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
    p=gscatter(x,y);
    hold on
    e=errorbar(x(1), y(1),Net(plotNet).random100.StdPath);
    e2=errorbar(x(1), y(1),Net(plotNet).random100.StdGlobalClust,'horizontal');
    errorbar(x(3), y(3),Net(plotNet).ordered100.StdPath);
    errorbar(x(3), y(3),Net(plotNet).ordered100.StdGlobalClust,'horizontal');
    % xlim([0.05 0.6])
    % ylim([2 16])
    text(x,y,{[num2str(Net(plotNet).sizeNetwork) 'node Random Nw'],'C. Elegans Nw', 'Human Nw',[num2str(Net(plotNet).sizeNetwork) 'node Ordered Nw'],'100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
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
% e=errorbar(logx(1), y(1),Net(plotNet).random100.StdPath);
% e2=errorbar(logx(1), y(1),Net(plotNet).random100.StdGlobalClust);
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

%Small World Prop
f2=figure;
for i = 1:plotNet
x=[Net(i).random100.AvgSmallWorldProp Net(i).ordered100.AvgSmallWorldProp cElegans.SmallWorldProp AgNW.SmallWorldProp];
p2=bar(x);
hold on
e=errorbar(x(1), Net(i).random100.StdSmallWorldProp);
e2=errorbar(x(2),Net(i).ordered100.StdSmallWorldProp);
% xlim([0.05 0.6])
% ylim([2 16])
xticklabels({[num2str(Net(i).sizeNetwork) 'node Random Nw'],[num2str(Net(i).sizeNetwork) 'node Ordered Nw'],'C. Elegans Nw','100nw','500nw','1000nw','2000nw'});
ylabel('Small World Prop');
hold on 
end 

fprintf('Figure 2 Complete \n');


f3=figure;
circuitRankRandom=[];
circuitRankOrdered=[];
for i = 1:plotNet
%Circuit Rank:
circuitRankRandom=[circuitRankRandom Net(i).random(1).CircuitRank];
circuitRankOrdered=[circuitRankOrdered Net(i).ordered(1).CircuitRank];
randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Random Nw'];
orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
end 

circuitRank=[circuitRankRandom circuitRankOrdered AgNW.CircuitRank];
p3=bar(circuitRank);
xticklabels([randomLabel orderedLabel {'100nw', '500nw', '1000nw','2000nw'}]);
ylabel('Circuit Rank');
hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles
for i=1:length(p2)  % iterate over number of bar objects
    hT=[hT text(p2(i).XData+p2(i).XOffset,p2(i).YData,num2str(p2(i).YData.'), ...
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
for i = 1:plotNet
%Circuit Rank:
PRandom=[PRandom Net(i).random100.AvgAvgPCoeff];
POrdered=[POrdered Net(i).ordered100.AvgAvgPCoeff];
MZRandom=[MZRandom Net(i).random100.AvgAvgMZ];
MZOrdered=[MZOrdered Net(i).ordered100.AvgAvgMZ];
randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Random Nw'];
orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
end 

PCoeff=[PRandom POrdered cElegans.avgP human.AvgP human.PLocalHubs human.PLocalHubs human.PConnectorHubs human.PConnectorHubs  mean(e100.Explore.GraphTheory.P), mean(e500.Explore.GraphTheory.P), mean(e1000.Explore.GraphTheory.P) mean(e2000.Explore.GraphTheory.P)];
MZ=[MZRandom MZOrdered cElegans.avgMZ human.AvgMZ human.MZHubs human.MZNonHubs human.MZHubs human.MZNonHubs  mean(e100.Explore.GraphTheory.MZ), mean(e500.Explore.GraphTheory.MZ), mean(e1000.Explore.GraphTheory.MZ) mean(e2000.Explore.GraphTheory.MZ)];

p4=gscatter(PCoeff,MZ);
%High PCoeff = Hubs / Central areas (Power et al., 2013)
text(PCoeff,MZ,{randomLabel{:}, orderedLabel{:}, 'C. Elegans Nw', 'Human Average', 'Human Connector Local Provincial Hub','Human Local Peripheral Node','Human Connector Hub','Human Satellite Connector', '100nw Avg', '500nw Avg', '1000nw Avg','2000nw Avg'},'NorthWest');
xlabel('Average Participant Coefficient Coefficient');
ylabel('Average Module z-Score');
p4.MarkerEdgeColor='b';
p4(:,1).MarkerEdgeColor='r';
p4.LineWidth=1.5;
hold on 

fprintf('Figure 4 Complete \n');

%Plot Guimera & Amaral rectangles:
for i=1:plotNet
f5(i)=figure;
guimera(network,Net,plotNet); %change network here
fprintf(['Figure 5 part ' num2str(i) ' Complete \n']);

end 

fprintf('Figure 5 Complete \n');

% %Communicability
% f6=figure;
% COMM=[Net(plotNet).random100.AvgCOMM(1) Net(plotNet).ordered100.AvgCOMM(1) AgNW.AvgCOMM];
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
BCRandom=[BCRandom Net(i).random100.AvgBC];
BCOrdered=[POrdered Net(i).ordered100.AvgBC];
BCRandomstd=[BCRandomstd Net(i).random100.StdBC];
BCOrderedstd=[BCOrderedstd Net(i).ordered100.StdBC];
randomLabel{i}=[num2str(Net(i).sizeNetwork) ' Random Nw'];
orderedLabel{i}=[num2str(Net(i).sizeNetwork) ' Ordered Nw'];
end 

BC=[BCRandom BCOrdered cElegans.avgBC AgNW.AvgBC];
stdBC=[BCRandomstd BCOrderedstd cElegans.stdBC AgNW.StdBC];
p7=bar(BC);
hold on
e=errorbar(BC, stdBC);
e.LineStyle='none';
% xlim([0.05 0.6])
% ylim([2 16])
xticklabels({[num2str(Net(i).sizeNetwork) 'node Random Nw'],[num2str(Net(i).sizeNetwork) 'node Ordered Nw'],'C. Elegans Nw','100nw','500nw','1000nw','2000nw'});
ylabel('Betweenness Centrality');
hold on 
 
 
fprintf('Figure 7 Complete \n');


%% SAVE GRAPHS 
print(f,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Clustering Coefficient vs Path Length all networks.pdf']);
print(f2,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Small World Prop all networks.pdf']);
print(f3,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Circuit Rank all networks.pdf']);
print(f4,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Participant Coefficient vs Module z-Score all networks.pdf']);
for i = 1:length(f5)
print(f5(i),'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Guimera PC vs Module z-Score ' num2str(Net(i).sizeNetwork) 'nw network.pdf']);
end 
% print(f6,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Communicability all networks.pdf']);
print(f7,'-painters','-dpdf','-bestfit','-r600',[fig_dir 'Betweenness Centrality all networks.pdf']);

end

%% FUNCTIONS
function guimera(network, Net,plotNet)
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
PCoeffRandom=Net(plotNet).random100.AvgPCoeff;
PCoeffOrdered= Net(plotNet).ordered100.AvgPCoeff;
MZRandom=Net(plotNet).random100.AvgMZ;
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
            r6=scatter(PCoeff(i),MZ(i),[255 204 153]./255,'filled');
            r6.MarkerEdgeColor='black';
        case 7
            r7=scatter(PCoeff(i),MZ(i),[0.7 0.7 0.7],'filled');
            r7.MarkerEdgeColor='black';
    end
end
ylim([-3 8]);
xlabel('Participant Coefficient')
ylabel('Within-Module Degree z-Score');
title([num2str(Net(plotNet).sizeNetwork) 'nw Network']);
hold on 
end





%% OLD CODE
%% Load three networks
% n100=load('Net_Sx_20_NoW100_03_25-2019_11_23_38_.mat');
% n500=load('Net_Sx_20_NoW500_03_30-2019_11_16_59_.mat');
% n1000=load('Net_Sx_20_NoW1000_06_06-2019_11_33_53_.mat');
%
% %Load three simulations:
% s100=load('Net_Sx_20_NoW100_03_25-2019_11_23_38_Zdenka_Constant_6SimsOnly_4_Sec_Vmax_1_20-May-2019.mat');
% s500=load('Net_Sx_20_NoW500_03_30-2019_11_16_59_Zdenka_Constant_1SimsOnly_4_Sec2Electrodes_Vmax_0.25_06-Jun-2019.mat');
% s1000=load('Net_Sx_20_NoW1000_06_06-2019_11_33_53_Zdenka_Constant_1SimsOnly_4_Sec2Electrodes_Vmax_0.25_06-Jun-2019.mat');
%
%
% %Set up 100nw network
% if exist('s100.SelSims','var') == 1
%       temp= s100.SelSims;%save simulations from test network if it hasn't been changed
% else
%       temp = n100.network.Simulations;%save simulations from test network if it has been changed
% end
%
% n100.network.Simulations=temp;
%
% %Set up 500nw Network
% if exist('s500.SelSims','var') == 1
%       temp= s500.SelSims;%save simulations from test network if it hasn't been changed
% else
%       temp = n500.network.Simulations;%save simulations from test network if it has been changed
% end
%
% n500.network.Simulations=temp;
%
% %Set up 1000nw Network
% if exist('s1000.SelSims','var') == 1
%       temp= s1000.SelSims;%save simulations from test network if it hasn't been changed
% else
%       temp = n1000.network.Simulations;%save simulations from test network if it has been changed
% end
%
% n1000.network.Simulations=temp;