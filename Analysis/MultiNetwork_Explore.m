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

explore_location='D:\alon_\Research\POSTGRAD\PhD\CODE\Data\Explore Analysis\';

%Load three explore analyses:
e100=load([explore_location 'Adrian_Net_Sx20_NoW100_0325-2019_112338__Sim_1_SourceElectrode_6_DrainElectrode_76_Exploration_Analysis_ Timestamp_400_06-Jun-2019.mat']);
e500=load([explore_location 'Adrian_Net_Sx20_NoW500_0330-2019_111659__Sim_1_SourceElectrode_18_DrainElectrode_492_Exploration_Analysis_ Timestamp_400_06-Jun-2019.mat']);
e1000=load([explore_location 'Adrian_Net_Sx20_NoW1000_0606-2019_113353__Sim_1_SourceElectrode_32_DrainElectrode_1000_Exploration_Analysis_ Timestamp_400_06-Jun-2019.mat']);


   
%% Analysis
x=[e100.Explore.GraphTheory.GlobalC1ust, e500.Explore.GraphTheory.GlobalC1ust, e1000.Explore.GraphTheory.GlobalC1ust];
y=[e100.Explore.GraphTheory.AvgPath, e500.Explore.GraphTheory.AvgPath, e1000.Explore.GraphTheory.AvgPath];
f=figure;
p=plot(x,y,'o-');
xlim([0.35 0.42])
text(x,y,{'100nw','500nw','1000nw'},'VerticalAlignment','bottom','HorizontalAlignment','left')
xlabel('Global Clustering Coefficient');
ylabel('Global Mean Path Length');
p.MarkerEdgeColor='r';
p.LineWidth=1.5;
