%Save Graph Representation
                    explore_location_origin='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Continuous DC\';

    e100=load([explore_location_origin 'Zdenka_Net_Sx20_NoW100_0930-2019_110329_Length_6_Disp_0_Sim_1_Source_99_Drain_36_Explore_Timestamp_200.mat']);
            e500=load([explore_location_origin 'Zdenka_Net_Sx20_NoW500_0330-2019_111659_Length_6_Disp_0_Sim_1_Source_498_Drain_435_Explore_Timestamp_200.mat']);
            e1000=load([explore_location_origin 'Zdenka_Net_Sx20_NoW1000_0606-2019_113353_Length_6_Disp_0_Sim_1_Source_997_Drain_408_Explore_Timestamp_200.mat']);
            e2000=load([explore_location_origin 'Zdenka_Net_Sx20_NoW2000_0618-2019_125103_Length_6_Disp_0_Sim_1_Source_1999_Drain_935_Explore_Timestamp_200.mat']);
     
G=cell(6,1);            
G{1}=e100.Explore.GraphView.Graph;
G{2}=e500.Explore.GraphView.Graph;
G{3}=e1000.Explore.GraphView.Graph;
G{4}=e2000.Explore.GraphView.Graph;

%% Load ANN
load('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Artificial Neural Network Graph Theory\500-Nodes\ANNGraphTheory.mat');

G{5}=ANN.Graph;


%% load C Elegans
load('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Organic Networks Connectomes\cElegansGraphTheory.mat');
G{6}=elegans.Graph;

%% load Watts-Strogatz
load('C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\Ordered+Random Data\Random_Ordered_Graphs_networks.mat');
G{7}=randomOrdered.Graph{1,1};%Beta = 0;
G{8}=randomOrdered.Graph{2,1};%Beta = 0;
G{9}=randomOrdered.Graph{3,1};%Beta = 0;
G{10}=randomOrdered.Graph{4,1};%Beta = 0;
G{11}=randomOrdered.Graph{1,21};%Beta = 1;
G{12}=randomOrdered.Graph{2,21};%Beta = 1;
G{13}=randomOrdered.Graph{3,21};%Beta = 1;
G{14}=randomOrdered.Graph{4,21};%Beta = 1;
G{15}=randomOrdered.Graph{1,11};%Beta = 0.5;
G{16}=randomOrdered.Graph{2,11};%Beta = 0.5;
G{17}=randomOrdered.Graph{3,11};%Beta = 0.5;
G{18}=randomOrdered.Graph{4,11};%Beta = 0.5;


type={'100nw','500nw','1000nw','2000nw','500NodeANN','cElegans','100WSB0','500WSB0','1000WSB0','2000WSB0','100WSB1','500WSB1','1000WSB1','2000WSB1','100WSB05','500WSB05','1000WSB05','2000WSB05'};

%%
save_directory='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\';
for i = 1:length(G)
    colors={'0049E9','00551A','AA0707','6A0055','FC90F9','FC9C03','0049E9','00551A','AA0707','6A0055','0049E9','00551A','AA0707','6A0055','0049E9','00551A','AA0707','6A0055'};
    colorsEdge={'4d85ff','00cc3d','f62323','cc00a3','fee6fe','fdc468','4d85ff','00cc3d','f62323','cc00a3','4d85ff','00cc3d','f62323','cc00a3','4d85ff','00cc3d','f62323','cc00a3'};
    clrsRGB{i}=hex2rgb(colors{i});
    clrsRGBEdge{i}=hex2rgb(colorsEdge{i});
fG(i)=figure;
pG(i)=plot(G{i});
pG(i).NodeLabel=[];
pG(i).NodeColor=clrsRGB{i};
pG(i).EdgeColor=clrsRGBEdge{i};
pG(i).MarkerSize=2;
pG(i).EdgeAlpha=0.2;
fG(i).CurrentAxes.Box='Off';
fG(i).CurrentAxes.XAxis.Color='none';
fG(i).CurrentAxes.YAxis.Color='none';
print(fG(i),'-painters','-dpng','-r600',[save_directory num2str(type{i}) '_Graphical_Representation.png']);
end 

% print(fG,'-painters','-dpdf','-bestfit','-r300',[save_directory num2str(network.Name) 'Simulation' num2str(simNum) '_SourceElectrode_' num2str(Explore.GraphView.ElectrodePosition(1)) '_DrainElectrode_' num2str(Explore.GraphView.ElectrodePosition(2)) '_ExploreGraph.pdf']);
        
   

