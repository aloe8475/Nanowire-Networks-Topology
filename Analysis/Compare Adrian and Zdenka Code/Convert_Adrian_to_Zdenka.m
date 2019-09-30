%% Convert Adrian Code to Zdenka Code (Connectivity):
clear all

computer=getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        loadPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Networks\Adrian Networks';        
        cd(loadPath)
        waitfor(msgbox('Select the Network saved data'));
        [FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
        count=1;
        f{count}=fullfile(PathName,FileName);
        load(f{count});
        SelSims=SelNet.Simulations;
        
    case '' %if on linux
        currentPath='/suphys/aloe8475/Documents/CODE/Analysis';
    case 'LAPTOP-S1BV3HR7'
        Sim='D:\alon_\Research\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Explore\100nw_03_25_2019 - size 20 length6 disp0\Net_Sx_20_NoW100_03_25-2019_11_23_38_Zdenka_Square_9SimsOnly_4_Sec_2Electrodes_Vmax_0.75_24-Jul-2019.mat';
        Net='D:\alon_\Research\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Net_Sx_20_NoW100_03_25-2019_11_23_38_.mat';
        %case '' %--- Add other computer paths (e.g. Mike)
end

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
temp=largestcomponent(adj_matrix);
adj_matrix=adj_matrix(temp,temp); %find largest connected component

%Number of Wires
number_of_wires=length(xa);
% number_of_junctions=LayoutSim.

%Theta:


%Connectivity Edge Position

% Connectivity Edge Position for Network without Simulation
% %combine Int:
% for i=1:length(LayoutNet.IndexInt)
%     for j =length(LayoutNet.IndexInt{i}):-1:1
%         if LayoutNet.IndexInt{i}(j)< i 
%             LayoutNet.IndexInt{i}(j)=[];
%             LayoutNet.XInt{i}(j)=[];
%             LayoutNet.YInt{i}(j)=[];
%         end 
%     end 
% end 
 
% index=[LayoutNet.IndexInt{:}];
% xi=[LayoutNet.XInt{:}];
% yi=[LayoutNet.YInt{:}];

% Connectivity Edge Position for Network with Simulation
xi=triu(LayoutSim.CX);
xi=xi(adj_matrix~=0);
xi=xi(xi~=0)';

yi=triu(LayoutSim.CY);
yi=yi(adj_matrix~=0);
yi=yi(yi~=0)';

length_x=15;
length_y=15;
%%
save(['AdriantoZdenka' num2str(number_of_wires) 'nw_simulation' num2str(currSim) '.mat']);

%%
% Connectivity.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
% Connectivity.filename=['AdriantoZdenka' num2str(number_of_wires) 'nw_simulation' num2str(currSim) '.mat'];
% Connectivity=getConnectivity(Connectivity);
% 
% Connectivity2.WhichMatrix       = 'nanoWires';    % 'nanoWires' \ 'randAdjMat'
% Connectivity2.filename='2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat';
% Connectivity2=getConnectivity(Connectivity2);
% 
% %% Snapshots 2 figure:
% 
% %Snapshots:
% % 
% % 
% for i =1:length(SelSims{currSim}.Time)
%     frame.Timestamp=SelSims{currSim}.Time;    
%     %Convert Current in wires to current in junctions:
%     currs=triu(SelSims{currSim}.Data.Currents{i}); % Find values of currents
%     abscurrs=abs(currs);
%     Imat=full(abscurrs); %full current matrix (instead of sparse double) + absolute value
%     Ilist=Imat(triu(adj_matrix)~=0); % current list for junctions that are 1 in the Adj matrix (only 1 directional)
%     frame.Current=Ilist;
%     %Convert Resistance in wires to resistance in junctions:
%     Rmat=triu(SelSims{currSim}.Data.Rmat{i});
%     Rlist=Rmat(triu(adj_matrix)~=0);
%     frame.Resistance = Rlist;
%     frame.Voltage=frame.Current.*frame.Resistance; %V=IR;
% %     frame.OnOrOff    = compPtr.comp.OnOrOff;
% %     frame.filamentState = compPtr.comp.filamentState;
%     snapshots{i}=frame;
% end

% SimulationOptions.ContactNodes=SelSims{currSim}.Electrodes.PosIndex;
% axesLimits.VoltageCbar=[0 1];%minimum voltage to maximum voltage
%     axesLimits.DissipationCbar = [0,5]; % (1pW*10^0 : 1pW*10^5)
%     axesLimits.CurrentArrowSacling = 10;
% whatToPlot = struct(...
%                         'Nanowires',    true, ...
%                         'Contacts',     true, ...
%                         'Dissipation',  true, ...
%                         'Currents',     false, ...
%                         'Voltages',     true  ...
%                         );
% snapshotFigure = snapshotToFigure(snapshots{floor(length(snapshots)/2)},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);  

% %% Video:
% %Snapshot options:
% TimeVector = SelSims{currSim}.Time';
% numIterations=length(TimeVector);
%     snapshotPeriod   = 4*1e-3; % (sec) make it a multiple integer of dt
%     snapshotStep     = ceil(snapshotPeriod /1e-3);
%     snapshotsIdx     = 1:snapshotStep:numIterations;
%     
% v = VideoWriter('networkMovie','Motion JPEG AVI');
% 
%         v.FrameRate = floor(1/snapshotPeriod/25);
%         
%         v.Quality = 100;
%         open(v);
% for i = 1 : length(snapshots)
%             progressBar(i,length(snapshots));
%             frameFig = snapshotToFigure(snapshots{i},SimulationOptions.ContactNodes,Connectivity,whatToPlot,axesLimits);
%             writeVideo(v,getframe(frameFig));
%             close(frameFig);
% end
%          close(v);