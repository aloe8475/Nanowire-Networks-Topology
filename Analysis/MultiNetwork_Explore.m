%Across Network Exploration:

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
computer=getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        explore_location='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\';
    case '' %if on linux
        explore_location='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/';
    case 'LAPTOP-S1BV3HR7'
        explore_location='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Explore Analysis\';
    %case '' %--- Add other computer paths (e.g. Mike)
end

%Load three explore analyses:
e100=load([explore_location 'Adrian_Net_Sx20_NoW100_0325-2019_112338__Sim_1_SourceElectrode_6_DrainElectrode_76_Exploration_Analysis_ Timestamp_400_18-Jun-2019.mat']);
e500=load([explore_location 'Adrian_Net_Sx20_NoW500_0330-2019_111659__Sim_1_SourceElectrode_18_DrainElectrode_492_Exploration_Analysis_ Timestamp_400_18-Jun-2019.mat']);
e1000=load([explore_location 'Adrian_Net_Sx20_NoW1000_0606-2019_113353__Sim_1_SourceElectrode_32_DrainElectrode_1000_Exploration_Analysis_ Timestamp_400_18-Jun-2019.mat']);

%% Random Graph Analysis:
%Create random graph with same number of vertices as 500nw network, and the
%avg degree of 500nw network.
A=createRandRegGraph(height(e500.Explore.GraphView.Nodes),round(full(mean(e500.Explore.GraphTheory.DEG))));
random.Adj=A;
ranGraph=graph(A);
random.Graph=ranGraph;
%Random Cluster
[random.Ci,random.Q] = community_louvain(A,1);
[random.GlobalClust,random.AvgLocalClust, random.Clust] = clustCoeff(A);
%Random Path Length
random.Path = path_length(A);
random.AvgPath=mean(random.Path);
%Circuit Rank
random.CircuitRank=numedges(ranGraph) - (numnodes(ranGraph) - 1);

%% Human Graph Analysis
human.GlobalClust=0.53;
human.AvgPath=2.49;
%Taken from (Achard et al., 2006)

%% Ordered Graph Analysis
n = 25;
B = delsq(numgrid('S',n));
G = graph(B,'omitselfloops');
ordered.Adj=adjacency(G);
%Random Cluster
[ordered.Ci,ordered.Q] = community_louvain(ordered.Adj,1);
[ordered.GlobalClust,ordered.AvgLocalClust, ordered.Clust] = clustCoeff(ordered.Adj);
%Random Path Length
ordered.Path = path_length(ordered.Adj);
ordered.AvgPath=mean(ordered.Path);
ordered.Graph=G;
%Circuit Rank
ordered.CircuitRank=numedges(G) - (numnodes(G) - 1);

%% AgNW 
%Circuit Rank
AgNW.CircuitRank=[e100.Explore.GraphTheory.CircuitRank e500.Explore.GraphTheory.CircuitRank e1000.Explore.GraphTheory.CircuitRank];
AgNW.GlobalClust=[e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath];
AgNW.AvgPath=[e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath];

%% Plot:
% Small World Analysis
x=[random.GlobalClust human.GlobalClust ordered.GlobalClust e100.Explore.GraphTheory.GlobalClust, e500.Explore.GraphTheory.GlobalClust, e1000.Explore.GraphTheory.GlobalClust ];
y=[random.AvgPath human.AvgPath ordered.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath];
f=figure;
p=gscatter(x,y);
% xlim([0.05 0.6])
% ylim([2 16])
text(x,y,{'500nw Random Nw','Human Nw','500nw Ordered Nw','100nw','500nw','1000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
xlabel('Global Clustering Coefficient');
ylabel('Global Mean Path Length');
p.MarkerEdgeColor='b';
p(:,1).MarkerEdgeColor='r';
p.LineWidth=1.5;

%Circuit Rank:
circuitRank=[random.CircuitRank human.CircuitRank ordered.CircuitRank AgNW.CircuitRank];