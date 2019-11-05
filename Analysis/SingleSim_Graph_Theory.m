%% Single Sim Graph Analysis
%-------------------------------------------------------------------------

% This script loads Graph Theory analyses performed in Network_Explore for
% 1 simulation of the same size network.

% Author: Alon Loeffler
% Date: 16/08/2019
% Version: 1.0

% Changelog:
% 22/07/2019 - Added Guimera and Amaral (2005) colored rectangular
% distribution of participant coefficient and module z score
% 26/06/2019 - Added Participant coefficient, 	 World Propensity and
% Betweenness centrality measures.
%-------------------------------------------------------------------------
dbstop if error

% addpath(genpath('../'));
close all;

set(0,'DefaultFigureVisible','on')


computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        dataPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Explore Analysis\';
    case '' %linux LVM
        dataPath='/import/silo2/aloe8475/Documents/CODE/Data/Explore Analysis/Time Delay Analysis';
        
    case 'LAPTOP-S1BV3HR7'
        dataPath='D:\alon_\Research\PhD\CODE\Data\Explore Analysis\';
end
cd(dataPath)
load_data_question=lower(input('Load data? N - None, D - Explore Data \n','s'));
if load_data_question=='d'
    clear -except currpath computer
    if ~exist('netCOMM','var')
        waitfor(msgbox('Select the Explore saved data'));
        [FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
        f=fullfile(PathName,FileName);
        load(f);
    end
end

Sim=currentSim;
currpath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis';
cd(currpath);
%% Initialise Variables:
fprintf('Running Analysis \n');

%if not thresholded:
%             Adj2=Explore.GraphView.AdjMat;
%if thresholded:
Adj2=Explore.GraphView.AdjMat(threshold,threshold);
%
%         Adj2=Explore.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
%         Adj2=Adj2(threshold,threshold);
%Find currents
currs=abs(Sim.Data.Currents{end});%Find the current at the end time
%                 if ~isempty(currs)
currs=currs(threshold,threshold);

[jj,ii,~]=find(tril(Adj2));
cc=zeros(1,length(jj));
for k=1:length(jj)
    cc(k)=currs(ii(k),jj(k));
end

% extract lower triangular part of Adjacency matrix of network
[jj,ii,~]=find(tril(Adj2));
cc2=zeros(1,length(jj));

%Find edges in Adj matrix that have current in them

for k=1:length(jj)
    cc2(k)=Explore.GraphTheory.networkThreshold(ii(k),jj(k));
end

% remove edges in adj matrix that don't have current
cc3=cc(logical(cc2));
meanCC3=mean(cc3);
stdCC3=std(cc3);
%                 end
%% COMMUNICABILITY
currs=abs(Sim.Data.Currents{end});%Find the current at the end time
%                 if ~isempty(currs)
currs=currs(threshold,threshold);


% Adj=Explore.GraphView.AdjMat(threshold,threshold);
COMM=Explore.GraphTheory.COMM(threshold,threshold);

%         if largestcomponent
[jj,ii,~]=find(tril(Adj2));%));
%         else
%             [jj,ii,~]=find(tril(Adj));
%         end
%
com=zeros(1,length(jj));

for k=1:length(jj)
    com(k)=COMM(ii(k),jj(k));
end

% extract lower triangular part of Adjacency matrix of network
% Adj2=Explore.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
% Adj2=Adj2(threshold,threshold);
[jj,ii,~]=find(tril(Adj2)); %extract lower triangle
com2=zeros(1,length(jj));

for k=1:length(jj)
    com2(k)=Explore.GraphTheory.networkThreshold(ii(k),jj(k)); %Graph.networkThreshold is just the same as the thresholded adj matrix.
end

%Find Graph.COMM in network
com3=com(logical(com2));
com3(com3==Inf)=0;
meanCom3=mean(com3);
stdCom3=std(com3);

com3(isnan(com3))=0;

%% CharPath, Distance & Global Efficiency
%     [Distance, GE] = efficiency_bin(Adj,0);
%     Distance(Distance==Inf)=0;
%     % Distance from Source:
%     DistSource=Distance(idx{1},:);
%     [CharPath,Efficiency,ecc,Radius,Diameter]=charpath(Distance);

%% Current at Nodes:

% Take the sum:
sourceElec=Explore.GraphView.ElectrodePosition(1);
drainElec=Explore.GraphView.ElectrodePosition(2);

sourceCurrent=max(sum(Explore.GraphView.currents(sourceElec,:)));
drainCurrent=min(sum(Explore.GraphView.currents(drainElec,:)));

idx=Explore.Electrodes.Source;
idx2=Explore.Electrodes.Drain;

%% Betweenness Centrality

BC=Explore.GraphTheory.BC(threshold);
sourceBC=BC(idx);
drainBC=BC(idx2);
%% Path Length


PathLength = path_length(Adj2);
PathLength(PathLength==Inf)=0;
AvgPath=mean(PathLength);
%
%% Degree

Degree=Explore.GraphTheory.DEG(threshold);
%             TestDegree=Degree;
sourceDEG=Degree(idx);
drainDEG=Degree(idx2);

%% Participation Coefficient
PCoeff=Explore.GraphTheory.P(threshold);
sourcePCoeff=PCoeff(idx);
drainPCoeff=PCoeff(idx2);


%% Module z-Score
MZ=Explore.GraphTheory.MZ(threshold);
sourceMZ=MZ(idx);
drainMZ=MZ(idx2);


%% Modularity
%             Module=Explore.GraphTheory.Modularity(threshold);
%             %NEED TO ADD MODULARITY

%% Clustering
Clust=Explore.GraphTheory.Clust(threshold);
sourceClust=Clust(idx);
drainClust=Clust(idx2);

fprintf('Generating Plots... \n');

%% PLOTTING
f=figure;
pTime=plot(Sim.Data.IDrain1);
hold on
xlabel('Timesteps (0.01s)')
ylabel('Drain Current (A)');
yyaxis right
pTime=plot(Sim.Data.VSource1);
ylabel('Source Voltage (V)');
    line([Explore.IndexTime Explore.IndexTime], get(gca,'YLim'),'Color','r','LineStyle','--');
    hold on;
%     clear cc3 com3 PCoeff MZ AvgPath Degree Clust BC



%% Plot Graph
pGraph = plot(Explore.GraphView.Graph);
highlight(pGraph,[Explore.Electrodes.Source,Explore.Electrodes.Drain],'NodeColor','g','MarkerSize',4);

%% Plot timeseries:
fTime=figure('Position',[0 0 1920 1080]);
plot(Sim.Time, Sim.Data.VSource1);
ylim([0 1.5]);

title([num2str(length(Explore.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim.SimInfo.MaxV) 'V | DC Continuous | Timeseries']);

%% Plot COMM Correlations at Edges:
fCOMMCurr=figure('Position',[0 0 1920 1080]);
count = 1;
    
    s=scatter(com3,cc3);
    xlabel('Communicability');

    ylabel('Current (A)');
                title([num2str(length(Explore.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | DC Continuous']);
     hold on
    [rCorrelation.COMM,pCorrelation.COMM]=corrcoef(com3,cc3);    