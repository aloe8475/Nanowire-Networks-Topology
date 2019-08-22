%% MultiSim Explore Analysis
%-------------------------------------------------------------------------

% This script loads Graph Theory analyses performed in Network_Explore_MultiSim_Same_Network for
% all simulations of the same size network.

% Author: Alon Loeffler
% Date: 16/08/2019
% Version: 1.0

% Changelog:
% 22/07/2019 - Added Guimera and Amaral (2005) colored rectangular
% distribution of participant coefficient and module z score
% 26/06/2019 - Added Participant coefficient, Small World Propensity and
% Betweenness centrality measures.
%-------------------------------------------------------------------------


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
load_data_question=lower(input('Load data? N - None, D - Explore Data \n','s'));
if load_data_question=='d'
    clear all;
    if ~exist('netCOMM','var')
        waitfor(msgbox('Select the Explore saved data'));
        [FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
        f=fullfile(PathName,FileName);
        load(f);
    end
end

%% NEED TO ADD ANALYSIS FOR DIFFERENT TIMES:

% For each time point
for j = 1:length(Explore)
    %For each Simulations
    thisExplore=Explore{j};
    thisThreshold=threshold{j};
    fprintf(num2str(j));
    for i = 1:length(thisExplore)
        if ~isempty(thisExplore{i})
            %% What we are doing here is finding the adj matrix, and finding the edges that have current flowing through them.
            %
            %if we want to extract the largest connected component:
            if j > 1
                largestcomponent=1;
            else
                largestcomponent=0;
            end
            
            if largestcomponent
                [bin,binsize] =conncomp(thisExplore{i}.GraphView.Graph);
                id = binsize(bin) == max(binsize);
                G = subgraph(thisExplore{i}.GraphView.Graph, id);
                thisExplore{i}.GraphView.NodeIndices(~id)=[];
                Adj=adjacency(G);
            else
                
                G=thisExplore{i}.GraphView.Graph;
                Adj=adjacency(G);
            end
            
            %
            %     %%
            %     %Find index of electrodes:
            [m, idx{i}]=max(thisExplore{i}.GraphView.NodeIndices==thisExplore{i}.GraphView.ElectrodePosition(1));
            [m, idx2{i}]=max(thisExplore{i}.GraphView.NodeIndices==thisExplore{i}.GraphView.ElectrodePosition(2));
            
            
            
            %connect source to drain on a copy of the adj matrix - THIS IS JOEL'S
            %CODE HE IS A LEGEND
            if j > 1
                toDelete = false(length(Adj),1);
                for count=1:length(Adj)
                    tempAdj=Adj;
                    g=graph(Adj);
                    sp=shortestpath(g,count,idx{i});
                    for count2=1:length(sp)-1
                        g=g.rmedge(sp(count2),sp(count2+1));
                    end
                    sp2=shortestpath(g,count,idx2{i});
                    if isempty(sp2)
                        Adj(count,:)=0;
                        Adj(:,count)=0;
                    end
                    clear tempAdj
                end
            end
            
            
            %     paths = allpaths(tempAdj, idx{i}, idx2{i});
            %
            % nodes = 1:length(Explore{i}.GraphView.NodeIndices);
            % for y = 1:size(paths,1)
            %     mpath = paths{y,1};
            %     mcost = paths{y,2}; %number of steps
            %     for q = 1:length(mcost)
            %         p = mpath(q,:);
            %         for n1=nodes
            %             for n2=p
            %                 if n1 == n2
            %                     nodes(nodes == n1) = [];
            %                     break
            %                 end
            %             end
            %         end
            %     end
            % end
            %
            % %remove nodes not in loop:
            % Adj(nodes,:)=0;
            % Adj(:,nodes)=0;
            
            %         first = true;
            %         for n = p
            %             if first
            %                 first = false;
            %             else
            %                 fprintf(' -> ');
            %             end
            %             fprintf('%d', n);
            %         end
            %         fprintf('  cost: %d\n',mcost(q));
            %         end
            % end
            %
            
            %% Conductance
            
            %% Currents
            Adj2{i}=Adj;
            %         Adj{i}=Explore{i}.GraphView.AdjMat(thisThreshold{i},thisThreshold{i});
            %
            %         Adj2{i}=Explore{i}.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
            %         Adj2{i}=Adj2{i}(thisThreshold{i},thisThreshold{i});
            %Find currents
            currs{i}=(abs(Sim{i}.Data.Currents{end}));
            
            currs{i}=currs{i}(thisThreshold{i},thisThreshold{i});
            
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
                cc2(k)=thisExplore{i}.GraphTheory.networkThreshold(ii(k),jj(k));
            end
            
            % remove edges in adj matrix that don't have current
            cc3{i}=cc(logical(cc2));
            meanCC3{i}=mean(cc3{i});
            stdCC3{i}=std(cc3{i});
            
            %% COMMUNICABILITY
            % Adj{i}=Explore{i}.GraphView.AdjMat(thisThreshold{i},thisThreshold{i});
            COMM{i}=thisExplore{i}.GraphTheory.COMM(thisThreshold{i},thisThreshold{i});
            %         if largestcomponent
            [jj,ii,~]=find(tril(Adj));%{i}));
            %         else
            %             [jj,ii,~]=find(tril(Adj{i}));
            %         end
            %
            com=zeros(1,length(jj));
            
            for k=1:length(jj)
                com(k)=COMM{i}(ii(k),jj(k));
            end
            
            % extract lower triangular part of Adjacency matrix of network
            % Adj2{i}=Explore{i}.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
            % Adj2{i}=Adj2{i}(thisThreshold{i},thisThreshold{i});
            [jj,ii,~]=find(tril(Adj2{i})); %extract lower triangle
            com2=zeros(1,length(jj));
            
            for k=1:length(jj)
                com2(k)=thisExplore{i}.GraphTheory.networkThreshold(ii(k),jj(k)); %Graph.networkThreshold is just the same as the thresholded adj matrix.
            end
            
            %Find Graph.COMM in network
            com3{i}=com(logical(com2));
            meanCom3{i}=mean(com3{i});
            stdCom3{i}=std(com3{i});
            
            %% CharPath, Distance & Global Efficiency
            %     [Distance{i}, GE{i}] = efficiency_bin(Adj,0);
            %     Distance{i}(Distance{i}==Inf)=0;
            %     % Distance from Source:
            %     DistSource{i}=Distance{i}(idx{1},:);
            %     [CharPath{i},Efficiency{i},ecc{i},Radius{i},Diameter{i}]=charpath(Distance{i});
            
            %% Current at Nodes:
            
            % Take the sum:
            sourceElec{i}=thisExplore{i}.GraphView.ElectrodePosition(1);
            drainElec{i}=thisExplore{i}.GraphView.ElectrodePosition(2);
            
            sourceCurrent{i}=max(sum(thisExplore{i}.GraphView.currents(sourceElec{i},:)));
            drainCurrent{i}=min(sum(thisExplore{i}.GraphView.currents(drainElec{i},:)));
            
            %% Betweenness Centrality
            BC{i}=thisExplore{i}.GraphTheory.BC(thisThreshold{i});
            sourceBC{i}=BC{i}(idx{i});
            drainBC{i}=BC{i}(idx2{i});
            
            %% Path Length
            PathLength{i} = path_length(Adj);
            PathLength{i}(PathLength{i}==Inf)=0;
            AvgPath{i}=mean(PathLength{i});
            
            %% Degree
            Degree{i}=thisExplore{i}.GraphTheory.DEG(thisThreshold{i});
            sourceDEG{i}=Degree{i}(idx{i});
            drainDEG{i}=Degree{i}(idx2{i});
            
            %% Participation Coefficient
            PCoeff{i}=thisExplore{i}.GraphTheory.P(thisThreshold{i});
            sourcePCoeff{i}=PCoeff{i}(idx{i});
            drainPCoeff{i}=PCoeff{i}(idx2{i});
            
            
            %% Module z-Score
            MZ{i}=thisExplore{i}.GraphTheory.MZ(thisThreshold{i});
            sourceMZ{i}=MZ{i}(idx{i});
            drainMZ{i}=MZ{i}(idx2{i});
            
            
            %% Clustering
            Clust{i}=thisExplore{i}.GraphTheory.Clust(thisThreshold{i});
            sourceClust{i}=Clust{i}(idx{i});
            drainClust{i}=Clust{i}(idx2{i});
            
            %         %% Modularity
            %         Mod{i}=thisExplore{i}.GraphTheory.Modularity(thisThreshold{i});
            %         sourceMod{i}=Mod{i}(idx{i});
            %         drainMod{i}=Mod{i}(idx2{i});
            %
            %
            
            %     if largestcomponent
            %     [jj,ii,~]=find(tril(Adj));%{i}));
            %     else
            %         [jj,ii,~]=find(tril(Adj{i}));
            %     end
            % %
            %     between=zeros(1,length(jj));
            %
            %     for k=1:length(jj)
            %         between(k)=BC{i}(ii(k),jj(k));
            %     end
            %
            %     % extract lower triangular part of Adjacency matrix of network
            %     % Adj2{i}=Explore{i}.GraphView.AdjMat;%convert 498x498 matrix to EdgeCData ~(1x6065)
            %     % Adj2{i}=Adj2{i}(thisThreshold{i},thisThreshold{i});
            %     [jj,ii,~]=find(tril(Adj2{i})); %extract lower triangle
            %     between2=zeros(1,length(jj));
            %
            %     for k=1:length(jj)
            %         between2(k)=Explore{i}.GraphTheory.networkThreshold(ii(k),jj(k)); %Graph.networkThreshold is just the same as the thresholded adj matrix.
            %     end
            %
            %     %Find Graph.COMM in network
            %     between3{i}=between(logical(between2));
            
        end
    end
    
    
    %%
    %Combine com3 and cc3 (network comm and network currents):
    netCOMM{j}=[com3{:}];
    netCurrs{j}=[cc3{:}];
    netsourceCurrent=abs([sourceCurrent{:}]);
    netdrainCurrent{j}=abs([drainCurrent{:}]);
    netsourceBC{j}=[sourceBC{:}];
    netdrainBC{j}=[drainBC{:}];
    netsourcePCoeff{j}=[sourcePCoeff{:}];
    netdrainPCoeff{j}=[drainPCoeff{:}];
    netsourceMZ{j}=[sourceMZ{:}];
    netdrainMZ{j}=[drainMZ{:}];
    netsourceDEG{j}=[sourceDEG{:}];
    netdrainDEG{j}=[drainDEG{:}];
    netSourceIdx{j}=[idx{:}];
    netDrainIdx{j}=[idx2{:}];
    netSourceElec{j}=[sourceElec{:}];
    netDrainElec{j}=[drainElec{:}];
    %     netDistSource{j}=[DistSource{:}];
    %Means and Stds
    meanCOMM{j}=[meanCom3{:}];
    stdCOMM{j}=[stdCom3{:}];
    meanCurrs{j}=[meanCC3{:}];
    stdCurrs{j}=[stdCC3{:}];
end

%% Plots:

%Plot Graph
% p = plot(Explore{end}.GraphView.Graph,'NodeLabel',Explore{end}.GraphView.NodeIndices);

%Plot Correlations at Edges:
f=figure('Position',[0 0 1920 1080]);
s=scatter(netCOMM,netCurrs);
h=lsline; %Linear Fit

%Polynomial Fits -------
hp=polyfit(netCOMM,netCurrs,2); %2nd Order Polynomial Fit
x2=min(netCOMM):0.25:max(netCOMM);
y2=polyval(hp,x2);
hold on
p2=plot(x2,y2,'g');

hp3=polyfit(netCOMM,netCurrs,3); %3rd Order Polynomial Fit
x3=min(netCOMM):0.25:max(netCOMM);
y3=polyval(hp3,x3);
hold on
p3=plot(x3,y3,'m');
% ------

h.Color='r';
xlabel('Communicability');
ylabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Random Electrode Placement']);
hold on
legend([p2,p3, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
[r.COMM,p.COMM]=corrcoef(netCOMM,netCurrs);

%% Plot correlation current vs Path Length

%Plot Correlations at Edges:
fcurr=figure('Position',[0 0 1920 1080]);
scurr=scatter(netDistSource,netCurrs);
hcurr=lsline; %Linear Fit

%Polynomial Fits -------
hpcurr=polyfit(netDistSource,netCurrs,2); %2nd Order Polynomial Fit
x2=min(netDistSource):1:max(netDistSource);
y2=polyval(hpcurr,x2);
hold on
p2curr=plot(x2,y2,'g');

hp3curr=polyfit(netDistSource,netCurrs,3); %3rd Order Polynomial Fit
x3=min(netDistSource):1:max(netDistSource);
y3=polyval(hp3curr,x3);
hold on
p3curr=plot(x3,y3,'m');
% ------

hcurr.Color='r';
xlabel('Path Length');
ylabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Random Electrode Placement']);
hold on
legend([p2curr,p3curr, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
[r.PATH,p.PATH]=corrcoef(netDistSource,netCurrs);

%% Plot correlation for mean and variance
fmean=figure('Position',[0 0 1920 1080]);
smean=scatter(meanCOMM,meanCurrs);
e=errorbar(meanCOMM,stdCOMM);
e.LineStyle='none';
e.Marker='o';
e2=errorbar(meanCurrs,stdCurrs);
e2.LineStyle='none';
e2.Marker='o';
hmean=lsline; %Linear Fit
% hpmean=polyfit(meanCOMM,netCurrs,3); %3rd Order Polynomial Fit
% x2=min(meanCOMM):0.25:max(meanCOMM);
% y2=polyval(hpmean,x2);
% hold on
% plot(x2,y2)
h.Color='r';
xlabel('Communicability');
ylabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Mean & STD | Random Electrode Placement']);
[r.COMM,p.COMM]=corrcoef(meanCOMM,meanCurrs);


%Plot Correlations at Source Nodes
% Betweenness Centrality
f1=figure('Position',[0 0 1920 1080]);
s1=scatter(netsourceBC,netsourceCurrent);
h1=lsline;
% hp=polyfit(netsourceBC,netsourceCurrent,1); %1st Order Polynomial Fit
h1.Color='r';
xlabel('Betweenness Centrality');
ylabel('Current (A) at Source');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Source Electrode']);
[r.BCsource]=corrcoef(netsourceBC,netsourceCurrent);
hold on

labelpoints(netsourceBC,netsourceCurrent,netSourceIdx);

% Degree
f2=figure('Position',[0 0 1920 1080]);
s2=scatter(netsourceDEG,netsourceCurrent);
h2=lsline;
h2.Color='r';
xlabel('Degree');
ylabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Source Electrode']);
[r.DEGsource]=corrcoef(netsourceDEG,netsourceCurrent);
hold on

labelpoints(netsourceDEG,netsourceCurrent,netSourceIdx);

% Participation Coefficient & Module Z
f3=figure('Position',[0 0 1920 1080]);
s3=scatter(netsourcePCoeff,netsourceCurrent);
h3=lsline;
h3.Color='r';
xlabel('Participant Coefficient');
ylabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Source Electrode']);
[r.Psource]=corrcoef(netsourcePCoeff,netsourceCurrent);
hold on
labelpoints(netsourcePCoeff,netsourceCurrent,netSourceIdx);



%Plot Correlations at Drain Nodes
% Betweenness Centrality
f1=figure('Position',[0 0 1920 1080]);
s=scatter(netdrainBC,netdrainCurrent);
h2=lsline;
h2.Color='r';
xlabel('Betweenness Centrality');
ylabel('Current (A) at Drain');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Drain Electrode']);
[r.BCdrain]=corrcoef(netdrainBC,netdrainCurrent);
hold on
labelpoints(netdrainBC,netdrainCurrent,netDrainIdx);


% Degree
f2=figure('Position',[0 0 1920 1080]);
s2=scatter(netdrainDEG,netdrainCurrent);
h2=lsline;
h2.Color='r';
xlabel('Degree');
ylabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Drain Electrode']);
[r.DEGdrain]=corrcoef(netdrainDEG,netdrainCurrent);
labelpoints(netdrainDEG,netdrainCurrent,netDrainIdx);

% Participation Coefficient & Module Z
% To do: scatter for each simulation independently and label drain from subgraph instead of
% drain from original graph
f3=figure('Position',[0 0 1920 1080]);
yyaxis left
s3=scatter(netdrainCurrent,netdrainPCoeff);
hold on
ConnectorHub=patch([0 0 2.2*1e-4 2.2*1e-4],[0.3 0.7 0.7 0.3],'blue','LineStyle','none','FaceAlpha',0.2);
xlim([0.2*1e-4 2.2*1e-4])
labelpoints(netdrainCurrent,netdrainPCoeff,netDrainIdx);
h3=lsline;
h3.Color='blue';
ylabel('Participant Coefficient')
yyaxis right
s3b=scatter(netdrainCurrent,netdrainMZ,'r');
labelpoints(netdrainCurrent,netdrainMZ,netDrainIdx);
h3b=lsline;
h3b.Color='r';
%Plot rectangles:


ylabel('Module z-Score')
xlabel('Current (A)');
title([num2str(length(Explore{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.Settings.Time) ' sec | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Drain Electrode']);
[r.Pdrain]=corrcoef(netdrainPCoeff,netdrainCurrent);

