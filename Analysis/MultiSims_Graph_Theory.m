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
    progressBar(j,length(Explore));
    for i = 1:length(thisExplore)
        
        % Store the 'class' of time for the current time.
        if ~isempty(Explore{j}{i})
            idxTime{j}.Time(i)=Explore{j}{i}.IndexTime;
            if j==length(Explore)
                class{j}{i}='Never';
            elseif j<=3
                class{j}{i}='Early';
            elseif j>3 & j <7
                class{j}{i}='Mid';
            elseif j>=7 & j<length(Explore)
                class{j}{i}='Late';
            end
            
        else
            idxTime{j}.Time(i)=NaN;
%             if j==length(Explore)
%                 if ~isnan(idxTime{j}.Time(i))
%                     class{j}{i}='Never';
%                 else
%                     class{j}{i}=class{j-1}{i};
%                 end
%             elseif j<=3
%                 if ~isnan(idxTime{j}.Time(i))
%                     class{j}{i}='Early';
%                 else
%                     class{j}{i}=class{j-1}{i};
%                 end
%             elseif j>3 & j <7
%                 if ~isnan(idxTime{j}.Time(i))
%                     class{j}{i}='Mid';
%                 else
%                     class{j}{i}=class{j-1}{i};
%                 end
%             elseif j>=7 & j<length(Explore)
%                 if ~isnan(idxTime{j}.Time(i))
%                     class{j}{i}='Late';
%                 else
%                     class{j}{i}=class{j-1}{i};
%                 end
%             end
            
        end
        
        if ~isempty(thisExplore{i})
            %% What we are doing here is finding the adj matrix, and finding the edges that have current flowing through them.
            %
            %if we want to extract the largest connected component:
            
            
            
            if j > 0
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
            
            
            
            %% connect source to drain on a copy of the adj matrix - THIS IS JOEL'S
            %CODE HE IS A LEGEND
            
            if j > 0
                toDelete = false(length(Adj),1);
                for count=1:length(Adj)
                    tempAdj=Adj;
                    g=graph(Adj);
                    sp=shortestpath(g,count,idx{i});
                    for count2=1:length(sp)-1
                        g=g.rmedge(sp(count2),sp(count2+1));
                    end
                    sp2=shortestpath(g,count,idx2{i});
                    if isempty(sp2) & height(g.Nodes)>2 %we do not want to break nodes with only one edge between them
                        Adj(count,:)=0;
                        Adj(:,count)=0;
                    end
                    clear tempAdj
                end
            end
            gRemovedEdges{j}{i}=graph(Adj);
            if largestcomponent
                [bin,binsize] =conncomp(gRemovedEdges{j}{i});
                id = binsize(bin) == max(binsize);
                gRemovedEdges{j}{i} = subgraph(gRemovedEdges{j}{i}, id);
                thisExplore{i}.GraphView.NodeIndices(~id)=[];
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
%             TestDegree{j}{i}=Degree{i};
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
    
    category{j}.NaN=sum(isnan(idxTime{j}.Time));
    category{j}.Never=sum(strcmp(class{j},'Never'));
    category{j}.Early=sum(strcmp(class{j},'Early'));
    category{j}.Mid=sum(strcmp(class{j},'Mid'));
    category{j}.Late=sum(strcmp(class{j},'Late'));
    %%
    %Combine com3 and cc3 (network comm and network currents):
    netDegree{j}=[Degree{:}];
    netPathLength{j}=[AvgPath{:}];
    netCOMM{j}=[com3{:}];
    netCurrs{j}=[cc3{:}];
    netPCoeff{j}=vertcat(PCoeff{:});
    netMZ{j}=vertcat(MZ{:});
    netClust{j}=vertcat(Clust{:});
    netBC{j}=[BC{:}];
%     netsourceCurrent{j}=abs([sourceCurrent{:}]);
%     netdrainCurrent{j}=abs([drainCurrent{:}]);
%     netsourceBC{j}=[sourceBC{:}];
%     netdrainBC{j}=[drainBC{:}];
%     netsourcePCoeff{j}=[sourcePCoeff{:}];
%     netdrainPCoeff{j}=[drainPCoeff{:}];
%     netsourceMZ{j}=[sourceMZ{:}];
%     netdrainMZ{j}=[drainMZ{:}];
%     netsourceDEG{j}=[sourceDEG{:}];
%     netdrainDEG{j}=[drainDEG{:}];
%     netSourceIdx{j}=[idx{:}];
%     netDrainIdx{j}=[idx2{:}];
%     netSourceElec{j}=[sourceElec{:}];
%     netDrainElec{j}=[drainElec{:}];
    %     netDistSource{j}=[DistSource{:}];
    %Means and Stds
    meanCOMM{j}=[meanCom3{:}];
    stdCOMM{j}=[stdCom3{:}];
    meanCurrs{j}=[meanCC3{:}];
    stdCurrs{j}=[stdCC3{:}];
    
    clear meanCom3 stdCom3 meanCC3 stdCC3 com3 cc3 sourceCurrent sourceClust sourceBC sourcePCoeff sourceMZ sourceDEG idx idx2 sourceElec drainBC drainClust drainCurrent drainDEG drainElec drainMZ drainPCoeff
end

% Firstly we compare different categories (early, mid, late & never) and
% then we do for each square pulse individually

%% Categorical Processing:
logicalCategory.Early=strcmp(class{end},'Early');
logicalCategory.Mid=strcmp(class{end},'Mid');
logicalCategory.Late=strcmp(class{end},'Late');
logicalCategory.Never=strcmp(class{end},'Never');

%Degree:
endTime.Degree{1}=[Degree{logicalCategory.Early}];
endTime.Degree{2}=[Degree{logicalCategory.Mid}];
endTime.Degree{3}=[Degree{logicalCategory.Late}];
endTime.Degree{4}=[Degree{logicalCategory.Never}];

fCategories=figure;
edges=[0:max([endTime.Degree{:}])/7:max([endTime.Degree{:}])];
for j = 1:length(endTime.Degree)
    hCat=histogram(endTime.Degree{j},edges);
    values(j,:)=hCat.Values;  
end
bCat=bar(values,'grouped');
xticklabels({'Early','Mid','Late','Never'});
ylabel('Frequency');
title('Degree (Categories)');

for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'DEG');
clear values Legend


%Path Length:
endTime.PathLength{1}=[AvgPath{logicalCategory.Early}];
endTime.PathLength{2}=[AvgPath{logicalCategory.Mid}];
endTime.PathLength{3}=[AvgPath{logicalCategory.Late}];
endTime.PathLength{4}=[AvgPath{logicalCategory.Never}];

fCat1=figure;
edges=[0:max([endTime.PathLength{:}])/7:max([endTime.PathLength{:}])];
for j = 1:length(endTime.PathLength)
    hCat1=histogram(endTime.PathLength{j},edges);
    values(j,:)=hCat1.Values;  
end
bCat1=bar(values,'grouped');
xticklabels({'Early','Mid','Late','Never'});
ylabel('Frequency');
title('Avg Path Length (Categories)');

for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'Avg Path');

clear values Legend

%% Clustering Coeff:
endTime.Clust{1}=vertcat(Clust{logicalCategory.Early});
endTime.Clust{2}=vertcat(Clust{logicalCategory.Mid});
endTime.Clust{3}=vertcat(Clust{logicalCategory.Late});
endTime.Clust{4}=vertcat(Clust{logicalCategory.Never});

fCat2=figure;
edges=[0:max(vertcat(endTime.Clust{:}))/7:max(vertcat(endTime.Clust{:}))];
for j = 1:length(endTime.Clust)
    hCat2=histogram(endTime.Clust{j},edges);
    values(j,:)=hCat2.Values;  
end
bCat2=bar(values,'grouped');
xticklabels({'Early','Mid','Late','Never'});
ylabel('Frequency');
title('Clustering Coefficient (Categories)');

for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'Clust');

clear values Legend


%netPCoeff
endTime.PCoeff{1}=vertcat(PCoeff{logicalCategory.Early});
endTime.PCoeff{2}=vertcat(PCoeff{logicalCategory.Mid});
endTime.PCoeff{3}=vertcat(PCoeff{logicalCategory.Late});
endTime.PCoeff{4}=vertcat(PCoeff{logicalCategory.Never});

fCat3=figure;
subplot(2,1,1)
edges=[0:max(vertcat(endTime.PCoeff{:}))/7:max(vertcat(endTime.PCoeff{:}))];
for j = 1:length(endTime.PCoeff)
    hCat3=histogram(endTime.PCoeff{j},edges);
    values(j,:)=hCat3.Values;  
end
bCat3=bar(values,'grouped');
xticklabels({'Early','Mid','Late','Never'});
ylabel('Frequency');
title('Participant Coefficient (Categories)');

for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'PCoeff');

clear values Legend

%netMZ
endTime.MZ{1}=vertcat(MZ{logicalCategory.Early});
endTime.MZ{2}=vertcat(MZ{logicalCategory.Mid});
endTime.MZ{3}=vertcat(MZ{logicalCategory.Late});
endTime.MZ{4}=vertcat(MZ{logicalCategory.Never});

subplot(2,1,2)
edges=[0:max(vertcat(endTime.MZ{:}))/7:max(vertcat(endTime.MZ{:}))];
for j = 1:length(endTime.MZ)
    hCat3=histogram(endTime.MZ{j},edges);
    values(j,:)=hCat3.Values;  
end
bCat3=bar(values,'grouped');
xticklabels({'Early','Mid','Late','Never'});
ylabel('Frequency');
title('Within-Module Degree z-Score (Categories)');

for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'MZ');



%COMM
fCat4=figure;
endTime.COMM{1}=vertcat(com3{logicalCategory.Early});
endTime.COMM{2}=vertcat(com3{logicalCategory.Mid});
endTime.COMM{3}=vertcat(com3{logicalCategory.Late});
endTime.COMM{4}=vertcat(com3{logicalCategory.Never});

edges=[0:max(vertcat(endTime.COMM{:}))/7:max(vertcat(endTime.COMM{:}))];
for j = 1:length(endTime.COMM)
    hCat4=histogram(endTime.COMM{j},edges);
    values(j,:)=hCat4.Values;  
end
bCat4=bar(values,'grouped');
xticklabels({'Early','Mid','Late','Never'});
ylabel('Frequency');
title('Communicability (Categories)');

for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'COMM');

clear values Legend


%% Timestamps:

%Scatter Colours:
clrs={'r','g','b','c','m','y','k',[0.75 0.2 0.25],[0.3 0.75 0.55], [0.6 0.2 0.1], [0.4 0.75 0.2]};

%Degree:
f1=figure;
edges=[0:max([netDegree{:}])/7:max([netDegree{:}])];
for j = 1:length(netCurrs)
    h1=histogram(netDegree{j},edges);
    values(j,:)=h1.Values;  
end
b1=bar(values,'grouped');
% b=plot(values,'o-')
for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'DEG');
xlabel('Time');
ylabel('Frequency');
title('Degree');
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});

clear values h1 Legend

%Path Length:
f2=figure;
edges=[0:max([netPathLength{:}])/7:max([netPathLength{:}])];
for j = 1:length(netCurrs)
    h1=histogram(netPathLength{j},edges);
    values(j,:)=h1.Values;
    
end
b2=bar(values,'grouped');
% b=plot(values,'o-')
for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'Avg Path');
xlabel('Time');
ylabel('Frequency');
title('Mean Path Length');
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});

clear values h1 Legend
%Clustering Coeff:
f3=figure;
edges=[0:max(vertcat(netClust{:}))/7:max(vertcat(netClust{:}))];
for j = 1:length(netCurrs)
    h1=histogram(netClust{j},edges);
    values(j,:)=h1.Values;
    
end
b3=bar(values,'grouped');
% b=plot(values,'o-')
for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'Clust');
xlabel('Time');
ylabel('Frequency');
title('Clustering Coefficient');
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});

clear values h1 Legend


%COMM
f4=figure;
edges=[0:max([netCOMM{:}])/7:max([netCOMM{:}])];
for j = 1:length(netCurrs)
    h1=histogram(netCOMM{j},edges);
    values(j,:)=h1.Values;
    
end
b4=bar(values,'grouped');
% b=plot(values,'o-')
for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'COMM');
xlabel('Time');
ylabel('Frequency');
title('Communicability');
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});

clear values h1 Legend

%netPCoeff
f5=figure;

subplot(2,1,1)
edges=[0:max(vertcat(netPCoeff{:}))/7:max(vertcat(netPCoeff{:}))];
for j = 1:length(netCurrs)
    h1=histogram(netPCoeff{j},edges);
    values(j,:)=h1.Values;
    
end
b5=bar(values,'grouped');
% b=plot(values,'o-')
for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'PCoeff');
xlabel('Time');
ylabel('Frequency');
title('Participant Coefficient');
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
clear values h1 Legend


%netMZ
subplot(2,1,2)
edges=[0:max(vertcat(netMZ{:}))/7:max(vertcat(netMZ{:}))];
for j = 1:length(netCurrs)
    h1=histogram(netMZ{j},edges);
    values(j,:)=h1.Values;
    
end
b6=bar(values,'grouped');
% b=plot(values,'o-')
for i =1:length(edges)
% set(b(i),'FaceColor',clrs{i})
   
    Legend{i}=strcat(num2str(edges(i)));
end
l=legend(Legend,'Location','NorthWest');
title(l,'MZ');
xlabel('Time');
ylabel('Frequency');
title('Within-Module Degree z-Score');
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
clear values h1 Legend                                                                                                                                     


%% Plots:

%Plot Graph
% p = plot(Explore{end}.GraphView.Graph,'NodeLabel',Explore{end}.GraphView.NodeIndices);


%% Plot NaN/Early/Mid/Late/Never:


%Histogram of Times sampled
% fnan=figure('Position',[0 0 1920 1080]);
% for i = 1:length(num)
% hist(idxTime{i}.Time)
% hold on
% end 



%Plot of Max Current (NaNs)
categoryMat=cell2mat(category);
plot([1:11],[categoryMat.NaN]./length(idxTime{1}.Time),'o-')
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
ylabel('% Reached Max Current');
xlabel('Square Pulse Time (mSec)'); 

% hold on 
% plot([1:11],[numMat.Never]./length(idxTime{1}.Time),'o-')
% plot([1:11],[numMat.Early]./length(idxTime{1}.Time),'o-')
% plot([1:11],[numMat.Mid]./length(idxTime{1}.Time),'o-')
% plot([1:11],[numMat.Late]./length(idxTime{1}.Time),'o-')

% % Plot timeseries:
figure;
plot(Sim{2}.Time, Sim{2}.Data.VSource1);
ylim([0 0.75]);
xlabel('Seconds')
ylabel('Source (V)');
title([num2str(length(Explore{2}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{2}.SimInfo.MaxV) 'V | Random Electrode Placement | Timeseries']);



%Plot Correlations at Edges:
f=figure('Position',[0 0 1920 1080]);
for i = 1:length(Explore)
    
    s(i)=scatter(netCOMM{i},netCurrs{i},[],clrs{i});
    % h=lsline; %Linear Fit
    
    % %Polynomial Fits -------
    % hp=polyfit(netCOMM{i},netCurrs{i},2); %2nd Order Polynomial Fit
    % x2=min(netCOMM{i}):0.25:max(netCOMM{i});
    % y2=polyval(hp,x2);
    % hold on
    % p2=plot(x2,y2,'g');
    %
    % hp3=polyfit(netCOMM{i},netCurrs{i},3); %3rd Order Polynomial Fit
    % x3=min(netCOMM{i}):0.25:max(netCOMM{i});
    % y3=polyval(hp3,x3);
    % hold on
    % p3=plot(x3,y3,'m');
    % % ------
    
    % h.Color='r';
    xlabel('Communicability');
    ylabel('Current (A)');
    title([num2str(length(Explore{i}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Random Electrode Placement']);
    hold on
    % legend([p2,p3, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
    [r{i}.COMM,p{i}.COMM]=corrcoef(netCOMM{i},netCurrs{i});
    
    Legend{i}=strcat([num2str(Explore{i}{1}.IndexTime) ' sec']);
end
legend(Legend)


% log10 Current:
flog=figure('Position',[0 0 1920 1080]);
for i = 1:length(Explore)
    slog(i)=scatter(netCOMM{i},log10(netCurrs{i}),[],clrs{i});
    % h=lsline; %Linear Fit
    
    % %Polynomial Fits -------
    % hp=polyfit(netCOMM{i},netCurrs{i},2); %2nd Order Polynomial Fit
    % x2=min(netCOMM{i}):0.25:max(netCOMM{i});
    % y2=polyval(hp,x2);
    % hold on
    % p2=plot(x2,y2,'g');
    %
    % hp3=polyfit(netCOMM{i},netCurrs{i},3); %3rd Order Polynomial Fit
    % x3=min(netCOMM{i}):0.25:max(netCOMM{i});
    % y3=polyval(hp3,x3);
    % hold on
    % p3=plot(x3,y3,'m');
    % % ------
    
    % h.Color='r';
    xlabel('Communicability');
    ylabel('log10 Current (A)');
    title([num2str(length(Explore{i}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Random Electrode Placement']);
    hold on
    % legend([p2,p3, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
    [r{i}.logCOMM,p{i}.logCOMM]=corrcoef(netCOMM{i},log10(netCurrs{i}));
    
    LegendLog{i}=strcat([num2str(Explore{i}{1}.IndexTime) ' sec']);
end
legend(LegendLog)

