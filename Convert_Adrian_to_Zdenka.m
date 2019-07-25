%% Convert Adrian Code to Zdenka Code (Connectivity):
load('D:\alon_\Research\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Explore\100nw_03_25_2019\Subgraph Communicability\Net_Sx_20_NoW100_03_25-2019_11_23_38_Zdenka_Square_9SimsOnly_4_Sec_2Electrodes_Vmax_0.75_24-Jul-2019.mat')
load('D:\alon_\Research\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Net_Sx_20_NoW100_03_25-2019_11_23_38_.mat')

currSim=input('Which Simulation do you want to convert? \n?');

LayoutSim=SelSims{currSim}.SelLayout;
LayoutNet=SelNet.LayOut;

%'xc', 'yc', 'xi', 'yi', 'xa', 'ya', 'xb', 'yb'
% xc=LayoutSim.CX;
% yc=LayoutSim.CY;
xa=full(diag(LayoutSim.X1))';
xb=full(diag(LayoutSim.X2))';
ya=full(diag(LayoutSim.Y1))';
yb=full(diag(LayoutSim.Y2))';

xc=(xa+xb)/2;
yc=(ya+yb)/2;
%Euclidean Distances:
wire_distances=pdist2(xc',yc');
%Adj Mat
adj_matrix=LayoutSim.AdjMat;

%Number of Wires
number_of_wires=length(xa);
% number_of_junctions=LayoutSim.

%Theta:


%Connectivity Edge Position
LayoutNet.XInt;
LayoutNet.YInt;
LayoutNet.IndexInt;

%combine Int:
for i=1:length(LayoutNet.IndexInt)
    for j =length(LayoutNet.IndexInt{i}):-1:1
        if LayoutNet.IndexInt{i}(j)< i 
            LayoutNet.IndexInt{i}(j)=[];
            LayoutNet.XInt{i}(j)=[];
            LayoutNet.YInt{i}(j)=[];
        end 
    end 
end 

index=[LayoutNet.IndexInt{:}];
xi=[LayoutNet.XInt{:}];
yi=[LayoutNet.YInt{:}];

length_x=20;
length_y=20;
%%
save('AdriantoZdenka100nw.mat');

%%
Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
Connectivity.filename='AdriantoZdenka100nw.mat';
Connectivity=getConnectivity(Connectivity);

Connectivity2.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
Connectivity2.filename='2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
Connectivity2=getConnectivity(Connectivity2);