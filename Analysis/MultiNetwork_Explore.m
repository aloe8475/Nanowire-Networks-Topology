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
e2000=load([explore_location 'Adrian_Net_Sx20_NoW2000_0618-2019_125103__Sim_1_SourceElectrode_158_DrainElectrode_1820_Exploration_Analysis_ Timestamp_400_19-Jun-2019.mat']);

%% Human Graph Analysis
human.GlobalClust=0.53;
human.AvgPath=2.49;
%Taken from (Achard et al., 2006)

%% Random and Ordered Graph Analysis
loadPath='D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis\Network Explore Functions\';
    if exist([loadPath 'Ordered_Random_Graphs_500nw.mat'], 'file')
        load([loadPath 'Ordered_Random_Graphs_500nw.mat']);
    else
        fprintf('Ordered and Random Graphs have not been created yet \n');
                fprintf('Creating New Graphs \n');
    [random, ordered]=createRandom_Ordered_Graphs(e500);
    end 
%% AgNW 
%Circuit Rank
AgNW.CircuitRank=[e100.Explore.GraphTheory.CircuitRank e500.Explore.GraphTheory.CircuitRank e1000.Explore.GraphTheory.CircuitRank e2000.Explore.GraphTheory.CircuitRank];
AgNW.GlobalClust=[e100.Explore.GraphTheory.GlobalClust, e500.Explore.GraphTheory.GlobalClust, e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
AgNW.AvgPath=[e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];

%% Plot:
% Small World Analysis
x=[random.GlobalClust human.GlobalClust ordered.GlobalClust e100.Explore.GraphTheory.GlobalClust, e500.Explore.GraphTheory.GlobalClust, e1000.Explore.GraphTheory.GlobalClust e2000.Explore.GraphTheory.GlobalClust];
y=[random.AvgPath human.AvgPath ordered.AvgPath e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath e2000.Explore.GraphTheory.AvgPath];
f=figure;
p=gscatter(x,y);
% xlim([0.05 0.6])
% ylim([2 16])
text(x,y,{'500node Random Nw','Human Nw','500node Ordered Nw','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
xlabel('Global Clustering Coefficient');
ylabel('Global Mean Path Length');
p.MarkerEdgeColor='b';
p(:,1).MarkerEdgeColor='r';
p.LineWidth=1.5;

%Small World Log
f1=figure;
logx=log10(x);
logy=log10(y);
p1=gscatter(x,logy);
% xlim([0.05 0.6])
% ylim([2 16])
text(x,logy,{'500node Random Nw','Human Nw','500node Ordered Nw','100nw','500nw','1000nw','2000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
xlabel('Global Clustering Coefficient');
ylabel('Log10 Global Mean Path Length');
p1.MarkerEdgeColor='b';
p1(:,1).MarkerEdgeColor='r';
p1.LineWidth=1.5;

%Circuit Rank:
circuitRank=[random.CircuitRank ordered.CircuitRank AgNW.CircuitRank];
f2=figure;
p2=bar(circuitRank);
xticklabels({'500node Random Nw', '500node Ordered Nw', '100nw', '500nw', '1000nw','2000nw'});
ylabel('Circuit Rank'); 

hAx=gca;            % get a variable for the current axes handle
hT=[];              % placeholder for text object handles
for i=1:length(p2)  % iterate over number of bar objects
  hT=[hT text(p2(i).XData+p2(i).XOffset,p2(i).YData,num2str(p2(i).YData.'), ...
                         'VerticalAlignment','bottom','horizontalalign','center')];
end