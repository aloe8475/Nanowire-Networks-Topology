%% IV Curve Plot
clear all;
close all;

%load data from MonolithicDemo.m
load('2019-04-04-1.45A_9.9secOff_20secTotal_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat')

%call snapshotIV function
[voltage, resistance, current,SimulationOptions]=snapshotIV(snapshots,Connectivity,SimulationOptions); %Calculate I-V Characteristics of the snapshots

%call other IV function
% To Do 
%% Snapshot I-V
 function [voltage, resistance, current,SimulationOptions] =snapshotIV(snapshots,Connectivity,SimulationOptions)
%Convert cell array to array
snapshotsarray=cell2mat(snapshots);

%save voltage, resistance and current at each edge at each snapshot in a struct
voltage.data=[snapshotsarray.Voltage]; % numEdges x length(snapshots)
resistance.data=[snapshotsarray.Resistance];
current.data= 1e6*(voltage.data./resistance.data);

%snapshots at the on/off nodes:
%find which edges are attached to the on/off nodes:
edgesIn=find(Connectivity.EdgeList(1,:)==SimulationOptions.ContactNodes(1)); %input node edges 1
edgesIn=[edgesIn; find(Connectivity.EdgeList(2,:)==SimulationOptions.ContactNodes(1))']; %input node edges 2
edgesOut=find(Connectivity.EdgeList(1,:)==SimulationOptions.ContactNodes(2))';%output node edges 1
edgesOut=[edgesOut; find(Connectivity.EdgeList(2,:)==SimulationOptions.ContactNodes(2))']; %output node edges 2

SimulationOptions.edgesIn=edgesIn;
SimulationOptions.edgesOut=edgesOut;

%calculate voltage, resistance and current at all the edges attached to the input node
voltage.in=voltage.data(edgesIn,:); 
resistance.in=resistance.data(edgesIn,:);
current.in=1e6*(voltage.in./resistance.in);

%calculate voltage, resistance and current at the edges attached to the output node
voltage.out=voltage.data(edgesOut,:);
resistance.out=resistance.data(edgesOut,:);
current.out=1e6*(voltage.out./resistance.out);

% %averages
% voltage.meanAcrossTime=mean(voltage.data,2); 
% voltage.meanAcrossedges=mean(voltage.data);
% current.meanAcrossTime=mean(current.data,2);
% current.meanAcrossedges=mean(current.data);

%% Plot I-V curve (Snapshots)
figure;
p1=plot(voltage.data,current.data); %this will plot voltage x current for each edge, for each time snapshot
ylabel('Current (nA)');
xlabel('Voltage (V)');
title('I-V Characteristics');

%plot I-V curve with input and output voltage&current
figure;
pOut=plot(voltage.out',current.out');
hold on
pIn=plot(voltage.in',current.in');
ylabel('Current (nA)');
xlabel('Voltage (V)');
title('I-V Curve at Input and Output edges ');
legend(['Input Edge: ' num2str(edgesIn(1))], ['Input Edge: ' num2str(edgesIn(2))],['Input Edge: ' num2str(edgesIn(3))], ...
    ['Output Edge: ' num2str(edgesOut(1))], ['Output Edge: ' num2str(edgesOut(2))], ['Output Edge: ' num2str(edgesOut(3))],...
    ['Output Edge: ' num2str(edgesOut(4))], ['Output Edge: ' num2str(edgesOut(5))]);

% %Plot mean across Time
% figure;
% p2=plot(voltage.meanAcrossTime,current.meanAcrossTime);
% ylabel('Current (nA)')
% xlabel('Voltage (V)');
% title('Mean I-V Characteristics Across Time');
% 
% %Plot mean across edges
% figure;
% p2=plot(voltage.meanAcrossedges,current.meanAcrossedges);
% ylabel('Current (nA)')
% xlabel('Voltage (V)');
% title('Mean I-V Characteristics Across edges');
% 
end 
