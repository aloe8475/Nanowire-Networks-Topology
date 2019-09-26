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

set(0,'DefaultFigureVisible','on')


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

%% Initialise Variables:
numSimulations=length(Explore{1});
gRemovedEdges=cell(length(Explore),numSimulations);
gOriginal=cell(length(Explore),numSimulations);
netCurrs=cell(length(Explore),numSimulations);
netDegree=cell(length(Explore),numSimulations);
netPathLength=cell(length(Explore),numSimulations);
netCOMM=cell(length(Explore),numSimulations);
netCurrs=cell(length(Explore),numSimulations);
netPCoeff=cell(length(Explore),numSimulations);
netMZ=cell(length(Explore),numSimulations);
netClust=cell(length(Explore),numSimulations);
netBC=cell(length(Explore),numSimulations);

%if Adrian's code
if ~exist('analysis_type','var')
    analysis_type='e';
    equalGroups=1;
end

% For each time point
for j = 1:length(Explore)
    
    fprintf(['\n ' num2str(j) '\n ']);
    %For each Simulation
    thisExplore=Explore{j};
    thisThreshold=threshold{j};
    for i = 1:length(thisExplore)
        
        % Store the 'class' of time for the current time.
        if ~isempty(Explore{j}{i})
            idxTime{j}.Time(i)=Explore{j}{i}.IndexTime;
            
            if analysis_type == 't'
                %BINS
                counter=[];
                if j==3
                    counter=1;
                    class{counter}{i}='Third Pulse';
                    classNums{counter}=j;
                    classTime{counter}=idxTime{j};
                elseif j==4
                    counter=2;
                    class{counter}{i}='Mid of Time Delay';
                    classNums{counter}=j;
                    classTime{counter}=idxTime{j};
                elseif j == 5
                    counter=3;
                    class{counter}{i}='End of Time Delay';
                    classNums{counter}=j;
                    classTime{counter}=idxTime{j};
                elseif j == length(Explore)
                    counter=4;
                    class{counter}{i}='Fourth Pulse';
                    classNums{counter}=j;
                    classTime{counter}=idxTime{j};
                end
            else
                if j<=(1/3)*length(Explore)
                    class{j}{i}='Early';
                    classNums{1}=[1:(1/3)*length(Explore)];
                elseif j>(1/3)*length(Explore) & j <(2/3)*length(Explore)
                    class{j}{i}='Mid';
                    classNums{2}=[round((1/3)*length(Explore)):round((2/3)*length(Explore))];
                elseif j>=(2/3)*length(Explore) & j<length(Explore)
                    class{j}{i}='Late';
                    classNums{3}=[round(2/3*length(Explore))+1:length(Explore)-1];
                elseif j == length(Explore)
                    class{j}{i}='Never';
                    classNums{4}=length(Explore);
                end
            end
            
        else
            if analysis_type == 'e'
                idxTime{j}.Time(i)=NaN;
                if j==length(Explore)
                    if ~isnan(idxTime{j}.Time(i))
                        class{j}{i}='Never';
                        classNums{4}=length(Explore);
                    else
                        class{j}{i}=class{j-1}{i};
                    end
                elseif j<=(1/3)*length(Explore)
                    if ~isnan(idxTime{j}.Time(i))
                        class{j}{i}='Early';
                        classNums{1}=[1:(1/3)*length(Explore)];
                    else
                        class{j}{i}=class{j-1}{i};
                    end
                elseif j>(1/3)*length(Explore) & j <(2/3)*length(Explore)
                    if ~isnan(idxTime{j}.Time(i))
                        class{j}{i}='Mid';
                        classNums{2}=[round((1/3)*length(Explore)):round((2/3)*length(Explore))];
                    else
                        class{j}{i}=class{j-1}{i};
                    end
                elseif j>=(2/3)*length(Explore) & j<length(Explore)
                    if ~isnan(idxTime{j}.Time(i))
                        class{j}{i}='Late';
                        classNums{3}=[round(2/3*length(Explore))+1:length(Explore)-1];
                    else
                        class{j}{i}=class{j-1}{i};
                    end
                end
            else
                idxTime{j}.Time(i)=NaN;
                counter=[];
                if j==3
                    counter = 1;
                    if ~isnan(idxTime{j}.Time(i))
                        class{counter}{i}='Third Pulse';
                    else
                        class{counter}{i}=class{counter-1}{i};
                    end
                elseif j==4
                    counter = 2;
                    if ~isnan(idxTime{j}.Time(i))
                        class{counter}{i}='Mid of Time Delay';
                    else
                        class{counter}{i}=class{counter-1}{i};
                    end
                elseif j==5
                    counter = 3;
                    if ~isnan(idxTime{j}.Time(i))
                        class{counter}{i}='End of Time Delay';
                    else
                        class{counter}{i}=class{counter-1}{i};
                    end
                elseif j==length(Explore)
                    counter = 4;
                    if ~isnan(idxTime{j}.Time(i))
                        class{counter}{i}='Fourth Pulse';
                    else
                        class{counter}{i}=class{counter-1}{i};
                    end
                end
            end
        end
        if ~isempty(thisExplore{i})
            
            %store graph:
            gOriginal{j,i}=thisExplore{i}.GraphView.Graph;
            
            %% What we are doing here is finding the adj matrix, and finding the edges that have current flowing through them.
            %
            %if we want to extract the largest connected component:
            largestcomponent=1;
            
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
            %     %Find index of electrodes: (as long as their
            if sum(thisExplore{i}.GraphView.NodeIndices==thisExplore{i}.GraphView.ElectrodePosition(1))>0
                [m, idx{i}]=max(thisExplore{i}.GraphView.NodeIndices==thisExplore{i}.GraphView.ElectrodePosition(1));
            else
                idx{i}=[];
            end
            if sum(thisExplore{i}.GraphView.NodeIndices==thisExplore{i}.GraphView.ElectrodePosition(2))>0
                [m, idx2{i}]=max(thisExplore{i}.GraphView.NodeIndices==thisExplore{i}.GraphView.ElectrodePosition(2));
            else
                idx2{i}=[];
            end
            
            
            %% connect source to drain on a copy of the adj matrix - THIS IS JOEL'S
            %CODE HE IS A LEGEND
            
            %             if j > 0
            if ~isempty(idx{i}) && ~isempty(idx2{i}) %if both electrodes aren't 'empty' we check if there is any extra paths on the main connected graph
                nonConnected(j,i)=0;
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
            else
                % MAKE LIST + DELETE LATER
                nonConnected(j,i)=1;
            end
            %             end
            % Save removed-edges Graphs
            gRemovedEdges{j,i}=graph(Adj);
            %Remove zero degree nodes:
            gRemovedEdges{j,i}=gRemovedEdges{j,i}.rmnode(find(gRemovedEdges{j,i}.degree==0));
            %             if largestcomponent
            %                 [bin,binsize] =conncomp(gRemovedEdges{j}{i});
            %                 id = binsize(bin) == max(binsize);
            %                 gRemovedEdges{j}{i} = subgraph(gRemovedEdges{j}{i}, id);
            %                 thisExplore{i}.GraphView.NodeIndices(~id)=[];
            %             end
            %
            
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
            if analysis_type=='t'
                if j>2
                    boolJ=true;
                else
                    boolJ=false;
                end
            else
                boolJ=true;
            end
            if boolJ %ONLY IF ANALYSIS_TYPE IS T
                currs{i}=(abs(Sim{i}.Data.Currents{idxTime{j}.Time(i)}));%Find the current at the currentTime
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
                com3{i}(com3{i}==Inf)=0;
                meanCom3{i}=mean(com3{i});
                stdCom3{i}=std(com3{i});
                
                com3{i}(isnan(com3{i}))=0;
                
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
                
                
                %% Modularity
                %             Module{i}=thisExplore{i}.GraphTheory.Modularity(thisThreshold{i});
                %             %NEED TO ADD MODULARITY
                
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
                
                %Graphs:
            end
            
        end %if ~isempty(thisExplore)
        
        %Combine com3 and cc3 (network comm and network currents):
        %     netMolularity{j}=[Module{:}];
        
        if analysis_type=='t'
            if j>2
                netDegree{j,i}=Degree{i};
                netPathLength{j,i}=AvgPath{i};
                netCOMM{j,i}=com3{i};
                netCOMM{j,i}(netCOMM{j,i}==Inf)=0; %NEED TO DISCUSS WITH MAC
                netCurrs{j,i}=cc3{i};
                netPCoeff{j,i}=PCoeff{i};
                netMZ{j,i}=MZ{i};
                netClust{j,i}=Clust{i};
                netBC{j,i}=BC{i};
            end
        else
            netDegree{j,i}=Degree{i};
            netPathLength{j,i}=AvgPath{i};
            netCOMM{j,i}=com3{i};
            netCOMM{j,i}(netCOMM{j,i}==Inf)=0; %NEED TO DISCUSS WITH MAC
            netCurrs{j,i}=cc3{i};
            netPCoeff{j,i}=PCoeff{i};
            netMZ{j,i}=MZ{i};
            netClust{j,i}=Clust{i};
            netBC{j,i}=BC{i};
        end
        progressBar(i,length(thisExplore));
        
    end %for loop
    
    if analysis_type == 'e'
        category{j}.NaN=sum(isnan(idxTime{j}.Time));
        category{j}.Never=sum(strcmp(class{j},'Never'));
        %     category{j}.NeverNotConnected=sum(nonConnected(j,:));
        %     category{j}.NeverConnected=category{j}.Never-category{j}.NeverNotConnected;
        category{j}.Early=sum(strcmp(class{j},'Early'));
        category{j}.Mid=sum(strcmp(class{j},'Mid'));
        category{j}.Late=sum(strcmp(class{j},'Late'));
    else
        counter =[];
        if j >2
            counter=j-2;
            category{counter}.NaN=sum(isnan(idxTime{j}.Time));
            %     category{j}.NeverNotConnected=sum(nonConnected(j,:));
            %     category{j}.NeverConnected=category{j}.Never-category{j}.NeverNotConnected;
            category{counter}.Third=sum(strcmp(class{counter},'Third Pulse'));
            category{counter}.MidPulse=sum(strcmp(class{counter},'Mid of Time Delay'));
            category{counter}.EndPulse=sum(strcmp(class{counter},'End of Time Delay'));
            category{counter}.Fourth=sum(strcmp(class{counter},'Fourth Pulse'));
        end
    end
end

%% Reshape Data
if equalGroups
% DEGREE
equal.Degree=[netDegree{:}];
equal.Degree=reshape(equal.Degree,[],3);
% COMM
equal.COMM=[netCOMM{:}];
equal.COMM=reshape(equal.COMM,[],3);
% PATH LENGTH
equal.PathLength=[netPathLength{:}];
equal.PathLength=reshape(equal.PathLength,[],3);
% PC & MZ
equal.PC=vertcat(netPCoeff{:});
equal.PC=reshape(equal.PC,[],3);
equal.MZ=vertcat(netMZ{:});
equal.MZ=reshape(equal.MZ,[],3);
% Clust
equal.Clust=vertcat(netClust{:});
equal.Clust=reshape(equal.Clust,[],3);
% CURRS
equal.Currs=[netCurrs{:}];
equal.Currs=reshape(equal.Currs,[],3);
end 

%% Unequal Data
% Firstly we compare different categories (early, mid, late & never) and
% then we do for each square pulse individually

%Non-Connected Graphs at each time:
NonConnectedTimes=sum(nonConnected,2);

%% Categorical Processing:
if analysis_type == 'e'
    logicalCategory.Early=strcmp(class{j},'Early');
    logicalCategory.Mid=strcmp(class{j},'Mid');
    logicalCategory.Late=strcmp(class{j},'Late');
    logicalCategory.Never=strcmp(class{j},'Never');
else
    logicalCategory.Third=strcmp(class{1},'Third Pulse');
    logicalCategory.MidTime=strcmp(class{2},'Mid of Time Delay');
    logicalCategory.EndTime=strcmp(class{3},'End of Time Delay');
    logicalCategory.Fourth=strcmp(class{4},'Fourth Pulse');
end
% logicalCategory.NeverAndNotConnected=nonConnected(j,:)==1;
% logicalCategory.NeverAndConnected=logical(logicalCategory.Never-logicalCategory.NeverAndNotConnected);

%% Graphs
if analysis_type == 'e'
    if equalGroups
        tempGraphs=reshape(gOriginal,[],3);
        tempRGraphs=reshape(gRemovedEdges,[],3);
        categories.originalGraphs.Early=tempGraphs(:,1);
        categories.originalGraphs.Mid=tempGraphs(:,2);
        categories.originalGraphs.Late=tempGraphs(:,3);
        categories.connectedGraphs.Early=tempRGraphs(:,1);
        categories.connectedGraphs.Mid=tempRGraphs(:,2);
        categories.connectedGraphs.Late=tempRGraphs(:,3);
    else
        for i=1:length(classNums) %times
            for j = 1:length(classNums{i})%number of simulations in each time
                if i ==1
                    categories.originalGraphs.Early{j}={gOriginal{classNums{i}(j),logicalCategory.Early}};
                    categories.connectedGraphs.Early{j}={gRemovedEdges{classNums{i}(j),logicalCategory.Early}};
                elseif i==2
                    categories.originalGraphs.Mid{j}={gOriginal{classNums{i}(j),logicalCategory.Mid}};
                    categories.connectedGraphs.Mid{j}={gRemovedEdges{classNums{i}(j),logicalCategory.Mid}};
                elseif i==3
                    categories.originalGraphs.Late{j}={gOriginal{classNums{i}(j),logicalCategory.Late}};
                    categories.connectedGraphs.Late{j}={gRemovedEdges{classNums{i}(j),logicalCategory.Late}};
                else
                    categories.originalGraphs.Never{j}={gOriginal{classNums{i}(j),logicalCategory.Never}};
                    categories.connectedGraphs.Never{j}={gRemovedEdges{classNums{i}(j),logicalCategory.Never}};
                end
            end
        end
        categories.originalGraphs.explanation="Each Cell Array represents the pulse number within the category";
        categories.connectedGraphs.explanation="Each Cell Array represents the pulse number within the category";
    end
else
    for i=1:length(classNums) %times
        for j = 1:length(classNums{i})%number of simulations in each time
            if i ==1
                categories.originalGraphs.Third{j}={gOriginal{i,logicalCategory.Third}};
                categories.connectedGraphs.Third{j}={gRemovedEdges{i,logicalCategory.Third}};
            elseif i==2
                categories.originalGraphs.MidTime{j}={gOriginal{i,logicalCategory.MidTime}};
                categories.connectedGraphs.MidTime{j}={gRemovedEdges{i,logicalCategory.MidTime}};
            elseif i==3
                categories.originalGraphs.EndTime{j}={gOriginal{i,logicalCategory.EndTime}};
                categories.connectedGraphs.EndTime{j}={gRemovedEdges{i,logicalCategory.EndTime}};
            else
                categories.originalGraphs.Fourth{j}={gOriginal{i,logicalCategory.Fourth}};
                categories.connectedGraphs.Fourth{j}={gRemovedEdges{i,logicalCategory.Fourth}};
            end
        end
    end
    categories.originalGraphs.explanation="Each Cell Array represents the pulse number";
    categories.connectedGraphs.explanation="Each Cell Array represents the pulse number";
end

%% Degree at endTime
if analysis_type == 'e'
    if equalGroups
        categories.Degree{1}=equal.Degree(:,1);
        categories.Degree{2}=equal.Degree(:,2);
        categories.Degree{3}=equal.Degree(:,3);
    else
        categories.Degree{1}=[Degree{logicalCategory.Early}];
        categories.Degree{2}=[Degree{logicalCategory.Mid}];
        categories.Degree{3}=[Degree{logicalCategory.Late}];
        categories.Degree{4}=[Degree{logicalCategory.Never}];
    end
else
    categories.Degree{1}=[netDegree{3,:}];
    categories.Degree{2}=[netDegree{4,:}];
    categories.Degree{3}=[netDegree{5,:}];
    categories.Degree{4}=[netDegree{6,:}];
end
%% Path Length:
if analysis_type == 'e'
    if equalGroups
        categories.PathLength{1}=equal.PathLength(:,1);
        categories.PathLength{2}=equal.PathLength(:,2);
        categories.PathLength{3}=equal.PathLength(:,3);
    else
        categories.PathLength{1}=[AvgPath{logicalCategory.Early}];
        categories.PathLength{2}=[AvgPath{logicalCategory.Mid}];
        categories.PathLength{3}=[AvgPath{logicalCategory.Late}];
        categories.PathLength{4}=[AvgPath{logicalCategory.Never}];
    end
else
    categories.PathLength{1}=[netPathLength{3,:}];
    categories.PathLength{2}=[netPathLength{4,:}];
    categories.PathLength{3}=[netPathLength{5,:}];
    categories.PathLength{4}=[netPathLength{6,:}];
end

%% Clustering Coeff:
if analysis_type == 'e'
    if equalGroups
        categories.Clust{1}=equal.Clust(:,1);
        categories.Clust{2}=equal.Clust(:,2);
        categories.Clust{3}=equal.Clust(:,3);
    else
        categories.Clust{1}=vertcat(Clust{logicalCategory.Early});
        categories.Clust{2}=vertcat(Clust{logicalCategory.Mid});
        categories.Clust{3}=vertcat(Clust{logicalCategory.Late});
        categories.Clust{4}=vertcat(Clust{logicalCategory.Never});
    end
else
    categories.Clust{1}=vertcat(netClust{3,:});
    categories.Clust{2}=vertcat(netClust{4,:});
    categories.Clust{3}=vertcat(netClust{5,:});
    categories.Clust{4}=vertcat(netClust{6,:});
end


%% netPCoeff
if analysis_type == 'e'
    if equalGroups
        categories.PCoeff{1}=equal.PC(:,1);
        categories.PCoeff{2}=equal.PC(:,2);
        categories.PCoeff{3}=equal.PC(:,3);
    else
        categories.PCoeff{1}=vertcat(PCoeff{logicalCategory.Early});
        categories.PCoeff{2}=vertcat(PCoeff{logicalCategory.Mid});
        categories.PCoeff{3}=vertcat(PCoeff{logicalCategory.Late});
        categories.PCoeff{4}=vertcat(PCoeff{logicalCategory.Never});
    end
else
    categories.PCoeff{1}=vertcat(netPCoeff{3,:});
    categories.PCoeff{2}=vertcat(netPCoeff{4,:});
    categories.PCoeff{3}=vertcat(netPCoeff{5,:});
    categories.PCoeff{4}=vertcat(netPCoeff{6,:});
end

%% netMZ
if analysis_type == 'e'
    if equalGroups
        categories.MZ{1}=equal.MZ(:,1);
        categories.MZ{2}=equal.MZ(:,2);
        categories.MZ{3}=equal.MZ(:,3);
    else
        categories.MZ{1}=vertcat(MZ{logicalCategory.Early});
        categories.MZ{2}=vertcat(MZ{logicalCategory.Mid});
        categories.MZ{3}=vertcat(MZ{logicalCategory.Late});
        categories.MZ{4}=vertcat(MZ{logicalCategory.Never});
    end
else
    categories.MZ{1}=vertcat(netMZ{3,:});
    categories.MZ{2}=vertcat(netMZ{4,:});
    categories.MZ{3}=vertcat(netMZ{5,:});
    categories.MZ{4}=vertcat(netMZ{6,:});
end

%% COMM
if analysis_type == 'e'
    if equalGroups
        categories.COMM{1}=equal.COMM(:,1);
        categories.COMM{2}=equal.COMM(:,2);
        categories.COMM{3}=equal.COMM(:,3);
    else
        % fCat4=figure('Position',[0 0 1920 1080]);
        categories.COMM{1}=[com3{logicalCategory.Early}];
        categories.COMM{2}=[com3{logicalCategory.Mid}];
        categories.COMM{3}=[com3{logicalCategory.Late}];
        categories.COMM{4}=[com3{logicalCategory.Never}];
    end
else
    categories.COMM{1}=[netCOMM{3,:}];
    categories.COMM{2}=[netCOMM{4,:}];
    categories.COMM{3}=[netCOMM{5,:}];
    categories.COMM{4}=[netCOMM{6,:}];
end

%% Currents - This needs to be slightly less than at the very end time point
% fCat5=figure;
if analysis_type == 'e'
       if equalGroups
        categories.Curr{1}=equal.Currs(:,1);
        categories.Curr{2}=equal.Currs(:,2);
        categories.Curr{3}=equal.Currs(:,3);
    else
    categories.Curr{1}=[cc3{logicalCategory.Early}];
    categories.Curr{2}=[cc3{logicalCategory.Mid}];
    categories.Curr{3}=[cc3{logicalCategory.Late}];
    categories.Curr{4}=[cc3{logicalCategory.Never}];
       end 
else
    categories.Curr{1}=[netCurrs{3,:}];
    categories.Curr{2}=[netCurrs{4,:}];
    categories.Curr{3}=[netCurrs{5,:}];
    categories.Curr{4}=[netCurrs{6,:}];
end
%
% edges=[1e-7:max([endTime.Curr{:}])/7:max([endTime.Curr{:}])];
% for j = 1:length(endTime.Curr)
%     hCat5=histogram(endTime.Curr{j},edges);
%     values(j,:)=hCat5.Values;
% end
% bCat5=bar(values,'grouped');
% xticklabels({'Early','Mid','Late/Never'});
% ylabel('Frequency');
% title('Current (Categories)');
%
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'Curr (A)');
%
% clear values Legend


%% Colors
%CHANGE COLOURS
lightblue=rgb('cyan');
blue=rgb('deep blue');
lightred=rgb('rose pink');
red=rgb('burnt red');

%Split colours into bins
rdif=[blue(1)-lightred(1)];
r=blue(1):(-rdif)/(length(categories.COMM)-1):lightred(1);
gdif=[blue(2)-lightred(2)];
g=blue(2):(-gdif)/(length(categories.COMM)-1):lightred(2);
bdif=[blue(3)-lightred(3)];
b=blue(3):(-bdif)/(length(categories.COMM)-1):lightred(3);

% %Scatter Colours:
clrs={[r; g; b]'};

%Split colours into bins
rdif=[lightblue(1)-blue(1)];
r2=lightblue(1):(-rdif)/(length(categories.Clust)-1):blue(1);
gdif=[lightblue(2)-blue(2)];
g2=lightblue(2):(-gdif)/(length(categories.Clust)-1):blue(2);
bdif=[lightblue(3)-blue(3)];
b2=lightblue(3):(-bdif)/(length(categories.Clust)-1):blue(3);

clrs2={[r2; g2; b2]'};
%% Plot Graph
% p = plot(Explore{end}.GraphView.Graph,'NodeLabel',Explore{end}.GraphView.NodeIndices);

%% Plot Subgraphs
% Example subgraphs for random Early, Mid and Late/Never:
% loop=1;
% while loop==1
%     randEarly=randi(length(Explore{2}));
%     count1=1;
%     randMid=randi(length(Explore{6}));
%     count2=1;
%     randLate=randi(length(Explore{10}));
%     count3=1;
%     if ~isempty(Explore{2}{randEarly}) && ~isempty(Explore{2}{randEarly}.GraphView.Nodes) && count1==1
%         fEarly=figure;
%         pEarly=plot(Explore{2}{randEarly}.GraphView.Graph,'NodeLabel',Explore{2}{randEarly}.GraphView.NodeIndices);
%         count1=2;
%         title('Early Subgraph');
%     end
%     if ~isempty(Explore{2}{randMid}) && ~isempty(Explore{2}{randMid}.GraphView.Nodes) && count2==1
%         fMid=figure;
%         pMid=plot(Explore{2}{randMid}.GraphView.Graph,'NodeLabel',Explore{2}{randMid}.GraphView.NodeIndices);
%         count2=2;
%         title('Mid Subgraph');
%     end
%     if ~isempty(Explore{2}{randLate}) && ~isempty(Explore{2}{randLate}.GraphView.Nodes) && count3==1
%         fLate=figure;
%         pLate=plot(Explore{2}{randLate}.GraphView.Graph,'NodeLabel',Explore{2}{randLate}.GraphView.NodeIndices);
%         count3=2;
%         title('Late/Never Subgraph');
%     end
%     if ~isempty(Explore{2}{randLate}) && ~isempty(Explore{2}{randLate}.GraphView.Nodes) && ~isempty(Explore{2}{randMid})  && ~isempty(Explore{2}{randMid}.GraphView.Nodes) && ~isempty(Explore{2}{randEarly})  && ~isempty(Explore{2}{randEarly}.GraphView.Nodes)
%         loop = loop+1;
%     end
% end

%% Plot of Max Current (NaNs)
if analysis_type == 'e'
    fnan=figure('Position',[0 0 1920 1080]);
    categoryMat=cell2mat(category);
    plot([1:length(Explore)],[categoryMat.NaN]./length(idxTime{1}.Time),'o-')
    xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
    ylabel('% Reached Max Current');
    xlabel('Square Pulse Time (mSec)');
else %Find the xticklabels
    %     fnan=figure('Position',[0 0 1920 1080]);
    % categoryMat=cell2mat(category);
    % plot([1:length(Explore)],[categoryMat.NaN]./length(idxTime{1}.Time),'o-')
    % xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
    % ylabel('% Reached Max Current');
    % xlabel('Square Pulse Time (mSec)');
end

% %% Plot timeseries:
% fTime=figure('Position',[0 0 1920 1080]);
% plot(Sim{3}.Time, Sim{3}.Data.VSource1);
% ylim([0 1.5]);
% xlabel('Seconds')
% ylabel('Source (V)');
% if analysis_type == 'e'
% title([num2str(length(Explore{3}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{3}.SimInfo.MaxV) 'V | Time+Categories | Timeseries']);
% else
% title([num2str(length(Explore{3}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{3}.SimInfo.MaxV) 'V | Time+Memory | Timeseries']);
% end
%% Plot COMM Correlations at Edges:
fCOMMCurr=figure('Position',[0 0 1920 1080]);
count = 1;
for i = 1:length(categories.COMM)
    s(i)=scatter(categories.COMM{i},categories.Curr{i},[],clrs{1}(i,:));
    %     h=lsline; %Linear Fit
    
    %     %Polynomial Fits -------
    %     hp=polyfit(categories.COMM{i},categories.Curr{i},2); %2nd Order Polynomial Fit
    %     x2=min(categories.COMM{i}):0.25:max(categories.COMM{i});
    %     y2=polyval(hp,x2);
    %     hold on
    %     p2=plot(x2,y2,'g');
    %
    %     hp3=polyfit(categories.COMM{i},categories.Curr{i},3); %3rd Order Polynomial Fit
    %     x3=min(categories.COMM{i}):0.25:max(categories.COMM{i});
    %     y3=polyval(hp3,x3);
    %     hold on
    %     p3=plot(x3,y3,'m');
    %     % % ------
    
    % h.Color='r';
    xlabel('Communicability');
    %     xlim([0 8])
    ylabel('Current (A)');
    for j = 1:5 %check the first five Explores just for naming purposes
        if ~isempty(Explore{i}{j})
            if analysis_type == 'e'
                title([num2str(length(Explore{i}{j}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | Time+Categories']);
            else
                title([num2str(length(Explore{i}{j}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | Time+Memory']);
                
            end
        end
    end
    hold on
    % legend([p2,p3, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
    [rCorrelation{i}.COMM,pCorrelation{i}.COMM]=corrcoef(categories.COMM{i},categories.Curr{i});
    %     Legend{i}=strcat([num2str(Explore{i}{3}.IndexTime) ' sec']);
end
% legend(Legend)
if analysis_type == 'e'
    legend({'Early','Mid','Late','Never'})
else
    legend({'Pulse 3','Mid Point of Time Delay','End Point of Time Delay', 'Pulse 4'})
end
%% log10 Current:
% flog=figure('Position',[0 0 1920 1080]);
% for i = 1:length(Explore)
%     slog(i)=scatter(netCOMM{i},log10(netCurrs{i}),[],clrs{i});
%     % h=lsline; %Linear Fit
%
%     % %Polynomial Fits -------
%     % hp=polyfit(netCOMM{i},netCurrs{i},2); %2nd Order Polynomial Fit
%     % x2=min(netCOMM{i}):0.25:max(netCOMM{i});
%     % y2=polyval(hp,x2);
%     % hold on
%     % p2=plot(x2,y2,'g');
%     %
%     % hp3=polyfit(netCOMM{i},netCurrs{i},3); %3rd Order Polynomial Fit
%     % x3=min(netCOMM{i}):0.25:max(netCOMM{i});
%     % y3=polyval(hp3,x3);
%     % hold on
%     % p3=plot(x3,y3,'m');
%     % % ------
%
%     % h.Color='r';
%     xlabel('Communicability');
%     ylabel('log10 Current (A)');
%     title([num2str(length(Explore{i}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{1}.SimInfo.MaxV) 'V | Random Electrode Placement']);
%     hold on
%     % legend([p2,p3, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
%     [r{i}.logCOMM,p{i}.logCOMM]=corrcoef(netCOMM{i},log10(netCurrs{i}));
%
%     LegendLog{i}=strcat([num2str(Explore{i}{1}.IndexTime) ' sec']);
% end
% legend(LegendLog)

%% Divide COMM into categories:
if analysis_type == 'e'
    fCOMMCat=figure('Position',[0 0 1920 1080]);
    scatterTemp1=0.9+0.2*rand(1,length(categories.COMM{1})); %create random noise
    scatterTemp2=1.9+0.2*rand(1,length(categories.COMM{2})); %create random noise
    scatterTemp3=2.9+0.2*rand(1,length(categories.COMM{3})); %create random noise
    
    scatterCat.Early=cell(1,length(categories.COMM{1}));
    scatterCat.Early(:)={'Early'};
    scatterCat.Mid=cell(1,length(categories.COMM{2}));
    scatterCat.Mid(:)={'Mid'};
    scatterCat.Late=cell(1,length(categories.COMM{3}));
    scatterCat.Late(:)={'Late'};
    if equalGroups
    gscatter([scatterTemp1 scatterTemp2 scatterTemp3]',vertcat(categories.COMM{1:3})',{scatterCat.Early{:},scatterCat.Mid{:},scatterCat.Late{:}}',[clrs2{1}],[],15,'off')
    title([num2str(length(Explore{3}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{3}.SimInfo.MaxV) 'V | Categories | Equal Groups']);
    else
        gscatter([scatterTemp1 scatterTemp2 scatterTemp3]',[categories.COMM{1:3}]',{scatterCat.Early{:},scatterCat.Mid{:},scatterCat.Late{:}}',[clrs2{1}],[],15,'off')
        title([num2str(length(Explore{3}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{3}.SimInfo.MaxV) 'V | Categories | Unequal Groups']);
    end 
    ylabel('Communicability');
    xticks([1 2 3 4]);
    xticklabels({'Early','Mid','Late'});
else
    fCOMMCat=figure('Position',[0 0 1920 1080]);
    scatterTemp1=0.9+0.2*rand(1,length(categories.COMM{1})); %create random noise
    scatterTemp2=1.9+0.2*rand(1,length(categories.COMM{2})); %create random noise
    scatterTemp3=2.9+0.2*rand(1,length(categories.COMM{3})); %create random noise
    scatterTemp4=3.9+0.2*rand(1,length(categories.COMM{4})); %create random noise
    
    scatterCat.Third=cell(1,length(categories.COMM{1}));
    scatterCat.Third(:)={'Pulse 3'};
    scatterCat.MidTime=cell(1,length(categories.COMM{2}));
    scatterCat.MidTime(:)={'Mid Point of Time Delay'};
    scatterCat.EndTime=cell(1,length(categories.COMM{3}));
    scatterCat.EndTime(:)={'End Point of Time Delay'};
    scatterCat.Fourth=cell(1,length(categories.COMM{4}));
    scatterCat.Fourth(:)={'Pulse 4'};
    gscatter([scatterTemp1 scatterTemp2 scatterTemp3 scatterTemp4]',cell2mat(categories.COMM)',{scatterCat.Third{:},scatterCat.MidTime{:},scatterCat.EndTime{:},scatterCat.Fourth{:}}',[clrs2{1}],[],15,'off')
    ylabel('Communicability');
    xticks([1 2 3 4 5]);
    xticklabels({'Pulse 3','Mid Point of Time Delay','End Point of Time Delay','Pulse 4'});
end


%% Clustering across categories:
% hist(categories.Clust{:})

%% Watson Strogatz Analysis:

if analysis_type=='e'
    clear cc;
    loadPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Watson Strogatz\';
    load([loadPath 'beta.mat'])
    watStr.beta=beta;
    load([loadPath 'cc.mat']);
    watStr.cc=cc;
    load([loadPath 'pl.mat']);
    watStr.pl=pl;
    
    
    fWatts=figure('Position',[0 0 1920 1080]);
    
    %Watss-Strogatz Colors:
    rdif=[red(1)-lightred(1)];
    r3=red(1):(-rdif)/100:lightred(1);
    gdif=[red(2)-lightred(2)];
    g3=red(2):(-gdif)/100:lightred(2);
    bdif=[red(3)-lightred(3)];
    b3=red(3):(-bdif)/100:lightred(3);
    
    s(1)=scatter(mean(watStr.cc),watStr.pl,[],[r3; g3; b3]');
    
    
    hold on
    for i=1:length(categories.Clust)
        %           s(i+1)=scatter(nanmean(categories.Clust{i}),categories.PathLength{i,j},[],clrs2{1}(i,:));
        s(i+1)=scatter(mean(categories.Clust{i}),nanmean(categories.PathLength{i}),[],clrs2{1}(i,:));
        hold on
        e=errorbar(mean(categories.Clust{i}),nanmean(categories.PathLength{i}),std(categories.Clust{i})/sqrt(length(categories.Clust{i})),'horizontal');
        e1=errorbar(mean(categories.Clust{i}),nanmean(categories.PathLength{i}),nanstd(categories.PathLength{i})/sqrt(length(categories.Clust{i})));
        e.Color=clrs2{1}(i,:);
        e1.Color=clrs2{1}(i,:);
    end
    l=legend(s,'Watss-Strogatz Distribution (Random to Ordered)','Early Subgraph','Mid Subgraph','Late Subgraph','location','NorthWest');
    ylabel('Avg Path Length');
    xlabel('Clustering Coefficient');
else
    %% Time Delay Clustering & Path Length
    clear cc;
    loadPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Watson Strogatz\';
    load([loadPath 'beta.mat'])
    watStr.beta=beta;
    load([loadPath 'cc.mat']);
    watStr.cc=cc;
    load([loadPath 'pl.mat']);
    watStr.pl=pl;
    
    
    fWatts=figure('Position',[0 0 1920 1080]);
    % Set up Legend
    %     Legend=cell(length(Explore),length(netClust));
    %
    %     s(1)=scatter(mean(watStr.cc),watStr.pl,[],[r3; g3; b3]');
    count=0;
    for i =3:length(Explore) %for each time window
        count=count+1;
        for j = 1:length(netClust) %for each simulation
            if ~isempty(netClust{i})
                s(count)=scatter3(nanmean(netClust{i,j}),netPathLength{i,j},j*0.05,[],clrs2{1}(count,:));
                hold on
            end
        end
    end
    axis vis3d
    ylabel('Avg Path Length');
    xlabel('Clustering Coefficient');
    zlabel('Time Delay');
    legend(s,{'Pulse 3','Mid Point of Time Delay','End Point of Time Delay','Pulse 4'});
    
    %% Time Delay Communicability
    count=0;
    catNums=[double(strcmp(class{1},'Third Pulse')); double(strcmp(class{2},'Mid of Time Delay'))*2; double(strcmp(class{3},'End of Time Delay'))*3; double(strcmp(class{4},'Fourth Pulse'))*4];
    for i =3:length(Explore) %for each time window
        count=count+1;
        for j = 1:length(netCOMM) %for each simulation
            if ~isempty(netCOMM{i})
                commTemp(count,j)=nanmean(netCOMM{i,j});
                %                 s=plot(0.05:0.05:5,commTemp);
                s(count)= scatter3(j*0.05,commTemp(count,j),catNums(count,j),[],clrs2{1}(count,:));
                hold on
            end
        end
    end
    
    
    %     for count = 1:length(s)
    %     s(count).LineWidth=1.5;
    %     s(count).Color=clrs2{1}(count,:);
    %     end
    axis vis3d
    xlabel('Time Delay');
    ylabel('Communicability');
    zticks([1 2 3 4])
    zticklabels({'Pulse 3','Mid Point of Time Delay','End Point of Time Delay','Pulse 4'});
end

%% 3D Surface Plots
x=commTemp;
y=catNums;
tri = delaunay(x,y);
z=[0.05:0.05:5;0.05:0.05:5;0.05:0.05:5;0.05:0.05:5];
h = trisurf(tri, x, y, z);
axis vis3d
colorbar EastOutside
zlabel('Time Delay (sec)');
yticks([1 2 3 4]);
xlabel('Communicability');
yticklabels({'Pulse 3','Mid Point of Time Delay','End Point of Time Delay','Pulse 4'});    %% ANOVAS
%% ANOVAS

if analysis_type=='e'
    %Make all the same length
    
    if ~equalGroups
    FigList = allchild(groot);
    
    maxlength = max(cellfun(@numel, categories.COMM));
    COMMtemp = cellfun(@(v) [v, nan(1, maxlength-numel(v))], categories.COMM, 'UniformOutput', false);
    else
        COMMtemp= categories.COMM;
    end 
    [ANOVA.COMM.p,ANOVA.COMM.AnovaTab,ANOVA.COMM.Stats] = anova1([COMMtemp{1}' COMMtemp{2}' COMMtemp{3}'],{'Early','Mid','Late'});
    [POSTHOC.COMM.c, POSTHOC.COMM.m, POSTHOC.COMM.h, POSTHOC.COMM.nms]=multcompare(ANOVA.COMM.Stats,'alpha',.05/3,'ctype','bonferroni');
    title('Communicability Post Hoc')
    
    FigHandle = setdiff(allchild(groot), FigList);
    FCOMMtable=FigHandle(1);
    FCOMMbox=FigHandle(2);
    
    clear FigHandle FigList
    
    %Path Length
    FigList = allchild(groot);
    
    maxlength = max(cellfun(@numel, categories.PathLength));
    PathLengthTemp = cellfun(@(v) [v, nan(1, maxlength-numel(v))], categories.PathLength, 'UniformOutput', false);
    
    [ANOVA.PathLength.p,ANOVA.PathLength.AnovaTab,ANOVA.PathLength.Stats] = anova1([PathLengthTemp{1}' PathLengthTemp{2}' PathLengthTemp{3}'],{'Early','Mid','Late'});
    [POSTHOC.PathLength.c, POSTHOC.PathLength.m, POSTHOC.PathLength.h, POSTHOC.PathLength.nms]=multcompare(ANOVA.PathLength.Stats,'alpha',.05/3,'ctype','bonferroni');
    title('Path Length Post Hoc')
    FigHandle = setdiff(allchild(groot), FigList);
    
    FPathtable=FigHandle(1);
    FPathbox=FigHandle(2);
    
    clear FigHandle FigList
    
    %Clustering Coeff
    FigList = allchild(groot);
    
    maxlength = max(cellfun(@numel, categories.Clust));
    ClustTemp = cellfun(@(v) [v, nan(1, maxlength-numel(v))], {categories.Clust{1}' categories.Clust{2}' categories.Clust{3}'}, 'UniformOutput', false);
    
    [ANOVA.Clust.p,ANOVA.Clust.AnovaTab,ANOVA.Clust.Stats] = anova1([ClustTemp{1}' ClustTemp{2}' ClustTemp{3}'],{'Early','Mid','Late'});
    title('Clustering Coeff ANOVA')
    [POSTHOC.Clust.c, POSTHOC.Clust.m, POSTHOC.Clust.h, POSTHOC.Clust.nms]=multcompare(ANOVA.Clust.Stats,'alpha',.05,'ctype','bonferroni');
    title('Clustering Coeff Post Hoc')
    FigHandle = setdiff(allchild(groot), FigList);
    
    FClusttable=FigHandle(1);
    FClustbox=FigHandle(2);
    
    clear FigHandle FigList
end
%% Save
save_directory='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Figures\Explore Analysis\Multiple Simulations - Same Network\';
MaxVoltage=num2str(Sim{1}.SimInfo.MaxV);
MaxVoltage=strrep(MaxVoltage,'.','');

% saveas(fCat,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_Degree_' date],'jpg');
% print(fCat,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_Degree_' date '.pdf']);
% saveas(fCat1,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_AvgPath_' date],'jpg');
% print(fCat1,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_AvgPath_' date '.pdf']);
% saveas(fCat2,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ClusteringCoeff_' date],'jpg');
% print(fCat2,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ClusteringCoeff_' date '.pdf']);
% saveas(fCat3,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_PCoeff_and_MZ_' date],'jpg');
% print(fCat3,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_PCoeff_and_MZ_' date '.pdf']);
% saveas(fCat4,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_Communicability_' date],'jpg');
% print(fCat4,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_Communicability_' date '.pdf']);
saveas(fnan,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_Time_%ReachedMaxCurrent_' date],'jpg');
print(fnan,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_Time_%ReachedMaxCurrent_' date '.pdf']);
saveas(fTime,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeSeries_' date],'jpg');
print(fTime,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeSeries_' date '.pdf']);
saveas(fCOMMCurr,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_Current vs Communicability Correlation_' date],'jpg');
print(fCOMMCurr,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TCurrent vs Communicability Correlation_' date '.pdf']);
saveas(FCOMMtable,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Table_COMM_' date],'jpg');
print(FCOMMtable,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Table_COMM_' date '.pdf']);
saveas(FCOMMbox,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Boxplot_COMM_' date],'jpg');
print(FCOMMbox,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Boxplot_COMM_' date '.pdf']);
saveas(FPathtable,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Table_PathLength_' date],'jpg');
print(FPathtable,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Table_PathLength_' date '.pdf']);
saveas(FPathbox,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Boxplot_PathLength_' date],'jpg');
print(FPathbox,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Boxplot_PathLength_' date '.pdf']);
saveas(FClusttable,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Table_ClustCoeff_' date],'jpg');
print(FClusttable,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Table_ClustCoeff_' date '.pdf']);
saveas(FClustbox,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_ANOVA_Boxplot_ClustCoeff_' date],'jpg');
print(FClustbox,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_PathLength vs ClustCoeff vs Watts-Strogatz' date '.pdf']);
saveas(fWatts,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_PathLength vs ClustCoeff vs Watts-Strogatz_' date],'jpg');
print(fWatts,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_PathLength vs ClustCoeff vs Watts-Strogatz_' date '.pdf']);
saveas(fCOMMCat,[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories_TimeCategories vs Communicability' date],'jpg');
print(fCOMMCat,'-painters','-dpdf','-bestfit','-r600',[save_directory num2str(length(Sim{1}.SelDomain)) 'nw_' num2str(length(network.Simulations)) 'sims_' MaxVoltage 'V_TimeCategories vs Communicability' date '.pdf']);


%% -----------------------------------------------------------------------
% %% Timestamps
%
% %Degree:
% f1=figure;
% edges=[0:max([netDegree{:}])/7:max([netDegree{:}])];
% for j = 1:length(netCurrs)
%     h1=histogram(netDegree{j},edges);
%     values(j,:)=h1.Values;
% end
% b1=bar(values,'grouped');
% % b=plot(values,'o-')
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'DEG');
% xlabel('Time');
% ylabel('Frequency');
% title('Degree');
% xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
%
% clear values h1 Legend
%
% %Path Length:
% f2=figure;
% edges=[0:max([netPathLength{:}])/7:max([netPathLength{:}])];
% for j = 1:length(netCurrs)
%     h1=histogram(netPathLength{j},edges);
%     values(j,:)=h1.Values;
%
% end
% b2=bar(values,'grouped');
% % b=plot(values,'o-')
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'Avg Path');
% xlabel('Time');
% ylabel('Frequency');
% title('Mean Path Length');
% xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
%
% clear values h1 Legend
% %Clustering Coeff:
% f3=figure;
% edges=[0:max(vertcat(netClust{:}))/7:max(vertcat(netClust{:}))];
% for j = 1:length(netCurrs)
%     h1=histogram(netClust{j},edges);
%     values(j,:)=h1.Values;
%
% end
% b3=bar(values,'grouped');
% % b=plot(values,'o-')
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'Clust');
% xlabel('Time');
% ylabel('Frequency');
% title('Clustering Coefficient');
% xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
%
% clear values h1 Legend
%
%
% %COMM
% f4=figure;
% edges=[0:max([netCOMM{:}])/7:max([netCOMM{:}])];
% for j = 1:length(netCurrs)
%     h1=histogram(netCOMM{j},edges);
%     values(j,:)=h1.Values;
%
% end
% b4=bar(values,'grouped');
% % b=plot(values,'o-')
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'COMM');
% xlabel('Time');
% ylabel('Frequency');
% title('Communicability');
% xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
%
% clear values h1 Legend
%
% %netPCoeff
% f5=figure;
%
% subplot(2,1,1)
% edges=[0:max(vertcat(netPCoeff{:}))/7:max(vertcat(netPCoeff{:}))];
% for j = 1:length(netCurrs)
%     h1=histogram(netPCoeff{j},edges);
%     values(j,:)=h1.Values;
%
% end
% b5=bar(values,'grouped');
% % b=plot(values,'o-')
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'PCoeff');
% xlabel('Time');
% ylabel('Frequency');
% title('Participant Coefficient');
% xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
% clear values h1 Legend
%
%
% %netMZ
% subplot(2,1,2)
% edges=[0:max(vertcat(netMZ{:}))/7:max(vertcat(netMZ{:}))];
% for j = 1:length(netCurrs)
%     h1=histogram(netMZ{j},edges);
%     values(j,:)=h1.Values;
%
% end
% b6=bar(values,'grouped');
% % b=plot(values,'o-')
% for i =1:length(edges)
% % set(b(i),'FaceColor',clrs{i})
%
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthWest');
% title(l,'MZ');
% xlabel('Time');
% ylabel('Frequency');
% title('Within-Module Degree z-Score');
% xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
% clear values h1 Legend


