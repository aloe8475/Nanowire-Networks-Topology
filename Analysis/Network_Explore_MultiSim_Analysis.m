%% MultiSim Explore Analysis
clear all;
close all;

computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        dataPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\';
    case ''
        dataPath='/suphys/aloe8475/Documents/CODE/Data/Explore Analysis/';
    case 'LAPTOP-S1BV3HR7'
        dataPath='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\';
end
cd(dataPath)

waitfor(msgbox('Select the Explore saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
f=fullfile(PathName,FileName);
load(f);


for i = 1:length(Explore)
    %% What we are doing here is finding the adj matrix, and finding the edges that have current flowing through them.
% %Find currents
% currs=currs(threshold,threshold); %currs(node_indices);%
% [j,i,~]=find(tril(Adj2));
% cc=zeros(1,length(j));
% for k=1:length(j)
%     cc(k)=currs(i(k),j(k));
% end
% 
% % extract lower triangular part of Adjacency matrix of network
% [j,i,~]=find(tril(Adj2));
% cc2=zeros(1,length(j));
% 
% %Find edges in Adj matrix that have current in them
% for k=1:length(j)
%     cc2(k)=Graph.networkThreshold(i(k),j(k));
% end
% % remove edges in adj matrix that don't have current
% cc3=cc(logical(cc2));
%     
    
    
    %% Currents
Adj{i}=Explore{i}.GraphView.AdjMat(threshold{i},threshold{i});
Adj2{i}=Explore{i}.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
Adj2{i}=Adj2{i}(threshold{i},threshold{i});
%Find currents
currs{i}=(abs(Sim{i}.Data.Currents{end}));

currs{i}=currs{i}(threshold{i},threshold{i});
[jj,ii,~]=find(tril(Adj2{i}));
cc=zeros(1,length(jj));
for k=1:length(jj)
    cc(k)=currs{i}(ii(k),jj(k));
end

% extract lower triangular part of Adjacency matrix of network
[jj,ii,~]=find(tril(Adj2{i}));
cc2=zeros(1,length(jj));

%Find edges in Adj matrix that have current in them

for k=1:length(jj)
    cc2(k)=Explore{i}.GraphTheory.networkThreshold(ii(k),jj(k));
end

% remove edges in adj matrix that don't have current
cc3{i}=cc(logical(cc2));
    

%% COMMUNICABILITY
% Adj{i}=Explore{i}.GraphView.AdjMat(threshold{i},threshold{i});
COMM{i}=Explore{i}.GraphTheory.COMM(threshold{i},threshold{i});

[jj,ii,~]=find(tril(Adj{i}));
com=zeros(1,length(j));

for k=1:length(jj)
    com(k)=COMM{i}(ii(k),jj(k));
end

% extract lower triangular part of Adjacency matrix of network
% Adj2{i}=Explore{i}.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
% Adj2{i}=Adj2{i}(threshold{i},threshold{i});
[jj,ii,~]=find(tril(Adj2{i})); %extract lower triangle
com2=zeros(1,length(jj)); 

for k=1:length(jj)
    com2(k)=Explore{i}.GraphTheory.networkThreshold(ii(k),jj(k)); %Graph.networkThreshold is just the same as the thresholded adj matrix.
end

%Find Graph.COMM in network
com3{i}=com(logical(com2));

end


%Combine com3 and cc3 (network comm and network currents):
netCOMM=[com3{:}];
netCurrs=[cc3{:}];

f=figure('Position',[0 0 1920 1080]);
f.PaperSize=[10,10];
s=scatter(netCOMM,netCurrs);
h1=lsline;
h1.Color='r';
xlabel('Communicability');
ylabel('Current (A)');
title('91nw | 50 Simulations | 1sec | 1V | Random Electrode Placement');
r=corrcoef(netCOMM,netCurrs);
% % For each simulation, calculate the correlation between communicability
% % and currents
% for i =1:length(Explore)
%     tmp=corrcoef(COMM{i},currs{i});
%     Rcoeff(i)=tmp(1,2);
% end 