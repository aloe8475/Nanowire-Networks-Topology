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
    
% For each time point
for j = 1:length(Explore)
    
    %For each Simulation
    thisExplore=Explore{j};
    thisThreshold=threshold{j};
    for i = 1:length(thisExplore)
        
        % Store the 'class' of time for the current time.
        if ~isempty(Explore{j}{i})
            idxTime{j}.Time(i)=Explore{j}{i}.IndexTime;
            if j<=(1/3)*length(Explore)
                class{j}{i}='Early';
            elseif j>(1/3)*length(Explore) & j <(2/3)*length(Explore)
                class{j}{i}='Mid';
            elseif j>=(2/3)*length(Explore) & j<length(Explore)
                class{j}{i}='Late';
            elseif j == length(Explore)
                class{j}{i}='Never';
            end
            
        else
            idxTime{j}.Time(i)=NaN;
            if j==length(Explore)
                if ~isnan(idxTime{j}.Time(i))
                    class{j}{i}='Never';
                else
                    class{j}{i}=class{j-1}{i};
                end
            elseif j<=(1/3)*length(Explore)
                if ~isnan(idxTime{j}.Time(i))
                    class{j}{i}='Early';
                else
                    class{j}{i}=class{j-1}{i};
                end
            elseif j>(1/3)*length(Explore) & j <(2/3)*length(Explore)
                if ~isnan(idxTime{j}.Time(i))
                    class{j}{i}='Mid';
                else
                    class{j}{i}=class{j-1}{i};
                end
            elseif j>=(2/3)*length(Explore) & j<length(Explore)
                if ~isnan(idxTime{j}.Time(i))
                    class{j}{i}='Late';
                else
                    class{j}{i}=class{j-1}{i};
                end        
            end 
        end
        
        if ~isempty(thisExplore{i})
            
            %store graph:
            gOriginal{j,i}=thisExplore{i}.GraphView.Graph;
            
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
            com3{i}(com3{i}==Inf)=0;
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
             
        end %if ~isempty(thisExplore)
        
            %Combine com3 and cc3 (network comm and network currents):
    %     netMolularity{j}=[Module{:}];
    netDegree{j,i}=Degree{i};
    netPathLength{j,i}=AvgPath{i};
    netCOMM{j,i}=com3{i};
    netCOMM{j,i}(netCOMM{j,i}==Inf)=0; %NEED TO DISCUSS WITH MAC
    netCurrs{j,i}=cc3{i};
    netPCoeff{j,i}=PCoeff{i};
    netMZ{j,i}=MZ{i};
    netClust{j,i}=Clust{i};
    netBC{j,i}=BC{i};
    end %for loop 
    
    
    category{j}.NaN=sum(isnan(idxTime{j}.Time));
    category{j}.Never=sum(strcmp(class{j},'Never'));
    category{j}.NeverNotConnected=sum(nonConnected(j,:));
    category{j}.NeverConnected=category{j}.Never-category{j}.NeverNotConnected;
    category{j}.Early=sum(strcmp(class{j},'Early'));
    category{j}.Mid=sum(strcmp(class{j},'Mid'));
    category{j}.Late=sum(strcmp(class{j},'Late'));
     
    progressBar(j,length(Explore));
end

% Firstly we compare different categories (early, mid, late & never) and
% then we do for each square pulse individually

%Non-Connected Graphs at each time:
NonConnectedTimes=sum(nonConnected,2);

%% Categorical Processing:
logicalCategory.Early=strcmp(class{end},'Early');
logicalCategory.Mid=strcmp(class{end},'Mid');
logicalCategory.Late=strcmp(class{end},'Late');
logicalCategory.Never=strcmp(class{end},'Never');
logicalCategory.NeverAndNotConnected=nonConnected(end,:)==1;
logicalCategory.NeverAndConnected=logical(logicalCategory.Never-logicalCategory.NeverAndNotConnected);

%% Graphs
categories.originalGraphs{1}={gOriginal{end,logicalCategory.Early}};
categories.originalGraphs{2}={gOriginal{end,logicalCategory.Mid}};
categories.originalGraphs{3}={gOriginal{end,logicalCategory.Late}};
categories.originalGraphs{4}={gOriginal{end,logicalCategory.NeverAndConnected}};

categories.connectedGraphs{1}={gRemovedEdges{end,logicalCategory.Early}};
categories.connectedGraphs{2}={gRemovedEdges{end,logicalCategory.Mid}};
categories.connectedGraphs{3}={gRemovedEdges{end,logicalCategory.Late}};
categories.connectedGraphs{4}={gRemovedEdges{end,logicalCategory.NeverAndConnected}};

%% Degree at endTime
categories.Degree{1}=[Degree{logicalCategory.Early}];
categories.Degree{2}=[Degree{logicalCategory.Mid}];
categories.Degree{3}=[Degree{logicalCategory.Late}];
% fCat=figure('Position',[0 0 1920 1080]);
% edges=[0:max([endTime.Degree{:}])/7:max([endTime.Degree{:}])];
% for j = 1:length(endTime.Degree)
%     hCat=histogram(endTime.Degree{j},edges);
%     values(j,:)=hCat.Values;
%     switch j
%         case 1
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 2
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 3
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 4
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%     end
% end
% bCat=bar(valuesPct,'grouped');
% xticklabels({['Early (' num2str(sum(values(1,:))) ')'],['Mid (' num2str(sum(values(2,:))) ')'],['Late (' num2str(sum(values(3,:))) ')'],['Never(' num2str(sum(values(1,:))) ')']});
% ylabel('Percentage');
% title('Degree (Categories)');

% for i =1:length(edges)
%     % set(b(i),'FaceColor',clrs{i})
%     
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthEast');
% title(l,'DEG');
% clear values Legend


%% Path Length:
categories.PathLength{1}=[AvgPath{logicalCategory.Early}];
categories.PathLength{2}=[AvgPath{logicalCategory.Mid}];
categories.PathLength{3}=[AvgPath{logicalCategory.Late}];

% fCat1=figure('Position',[0 0 1920 1080]);
% edges=[0:max([endTime.PathLength{:}])/7:max([endTime.PathLength{:}])];
% for j = 1:length(endTime.PathLength)
%     hCat1=histogram(endTime.PathLength{j},edges);
%     values(j,:)=hCat1.Values;
%     switch j
%         case 1
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 2
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 3
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 4
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%     end
% end
% bCat1=bar(valuesPct,'grouped');
% xticklabels({['Early (' num2str(sum(values(1,:))) ')'],['Mid (' num2str(sum(values(2,:))) ')'],['Late (' num2str(sum(values(3,:))) ')'],['Never (' num2str(sum(values(1,:))) ')']});
% ylabel('Percentage');
% title('Avg Path Length (Categories)');
% 
% for i =1:length(edges)
%     % set(b(i),'FaceColor',clrs{i})
%     
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthEast');
% title(l,'Avg Path');
% 
% clear values Legend

%% Clustering Coeff:
categories.Clust{1}=vertcat(Clust{logicalCategory.Early});
categories.Clust{2}=vertcat(Clust{logicalCategory.Mid});
categories.Clust{3}=vertcat(Clust{logicalCategory.Late});

% fCat2=figure('Position',[0 0 1920 1080]);
% edges=[0:max(vertcat(endTime.Clust{:}))/7:max(vertcat(endTime.Clust{:}))];
% for j = 1:length(endTime.Clust)
%     hCat2=histogram(endTime.Clust{j},edges);
%     values(j,:)=hCat2.Values;
%     switch j
%         case 1
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 2
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 3
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 4
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%     end
% end
% bCat2=bar(valuesPct,'grouped');
% xticklabels({['Early (' num2str(sum(values(1,:))) ')'],['Mid (' num2str(sum(values(2,:))) ')'],['Late (' num2str(sum(values(3,:))) ')'],['Never (' num2str(sum(values(1,:))) ')']});
% ylabel('Percentage');
% title('Clustering Coefficient (Categories)');
% 
% for i =1:length(edges)
%     % set(b(i),'FaceColor',clrs{i})
%     
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthEast');
% title(l,'Clust');
% 
% clear values Legend


%% netPCoeff
categories.PCoeff{1}=vertcat(PCoeff{logicalCategory.Early});
categories.PCoeff{2}=vertcat(PCoeff{logicalCategory.Mid});
categories.PCoeff{3}=vertcat(PCoeff{logicalCategory.Late});

% fCat3=figure('Position',[0 0 1920 1080]);
% subplot(2,1,1)
% edges=[0:max(vertcat(endTime.PCoeff{:}))/7:max(vertcat(endTime.PCoeff{:}))];
% for j = 1:length(endTime.PCoeff)
%     hCat3=histogram(endTime.PCoeff{j},edges);
%     values(j,:)=hCat3.Values;
%     switch j
%         case 1
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 2
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 3
%             valuesPctvaluesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 4
%             values(j,:) = (values(j,:)./sum(values(j,:)))*100;
%     end
% end
% bCat3=bar(valuesPct,'grouped');
% xticklabels({['Early (' num2str(sum(values(1,:))) ')'],['Mid (' num2str(sum(values(2,:))) ')'],['Late (' num2str(sum(values(3,:))) ')'],['Never (' num2str(sum(values(1,:))) ')']});
% ylabel('Percentage');
% title('Participant Coefficient (Categories)');
% 
% for i =1:length(edges)
%     % set(b(i),'FaceColor',clrs{i})
%     
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthEast');
% title(l,'PCoeff');
% 
% clear values Legend

%% netMZ
categories.MZ{1}=vertcat(MZ{logicalCategory.Early});
categories.MZ{2}=vertcat(MZ{logicalCategory.Mid});
categories.MZ{3}=vertcat(MZ{logicalCategory.Late});

% subplot(2,1,2)
% edges=[0:max(vertcat(endTime.MZ{:}))/7:max(vertcat(endTime.MZ{:}))];
% for j = 1:length(endTime.MZ)
%     hCat3=histogram(endTime.MZ{j},edges);
%     values(j,:)=hCat3.Values;
%     switch j
%         case 1
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 2
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 3
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 4
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%     end
% end
% bCat3a=bar(valuesPct,'grouped');
% xticklabels({['Early (' num2str(sum(values(1,:))) ')'],['Mid (' num2str(sum(values(2,:))) ')'],['Late (' num2str(sum(values(3,:))) ')'],['Never (' num2str(sum(values(1,:))) ')']});
% ylabel('Percentage');
% title('Within-Module Degree z-Score (Categories)');
% 
% for i =1:length(edges)
%     % set(b(i),'FaceColor',clrs{i})
%     
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthEast');
% title(l,'MZ');
% 
% clear values Legend

%% COMM
% fCat4=figure('Position',[0 0 1920 1080]);
categories.COMM{1}=[com3{logicalCategory.Early}];
categories.COMM{2}=[com3{logicalCategory.Mid}];
categories.COMM{3}=[com3{logicalCategory.Late}];
% for j = 1:length(endTime.COMM)
%     endTime.COMM{j}(endTime.COMM{j}==Inf)=0;
%     edges=[0:max([endTime.COMM{:}])/7:max([endTime.COMM{:}])];
%     hCat4=histogram(endTime.COMM{j},edges);
%     values(j,:)=hCat4.Values;
%     switch j
%         case 1
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 2
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 3
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%         case 4
%             valuesPct(j,:) = (values(j,:)./sum(values(j,:)))*100;
%     end
% end
% 
% bCat4=bar(valuesPct,'grouped');
% xticklabels({['Early (' num2str(sum(values(1,:))) ')'],['Mid (' num2str(sum(values(2,:))) ')'],['Late (' num2str(sum(values(3,:))) ')'],['Never (' num2str(sum(values(1,:))) ')']});
% ylabel('Percentage');
% title('Communicability (Categories)');
% 
% for i =1:length(edges)
%     % set(b(i),'FaceColor',clrs{i})
%     
%     Legend{i}=strcat(num2str(edges(i)));
% end
% l=legend(Legend,'Location','NorthEast');
% title(l,'COMM');
% 
% clear values Legend

%Currents - This needs to be slightly less than at the very end time point
% fCat5=figure;
categories.Curr{1}=[cc3{logicalCategory.Early}];
categories.Curr{2}=[cc3{logicalCategory.Mid}];
categories.Curr{3}=[cc3{logicalCategory.Late}];
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
lightblue=rgb('Light Blue');
blue=rgb('Blue');
lightred=rgb('rose pink');
red=rgb('burnt red');

%Split colours into bins
rdif=[lightblue(1)-blue(1)];
r=lightblue(1):(-rdif)/(length(categories.COMM)-1):blue(1);
gdif=[lightblue(2)-blue(2)];
g=lightblue(2):(-gdif)/(length(categories.COMM)-1):blue(2);
bdif=[lightblue(3)-blue(3)];
b=lightblue(3):(-bdif)/(length(categories.COMM)-1):blue(3);

% %Scatter Colours:
clrs={[r; g; b]'};
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
fnan=figure('Position',[0 0 1920 1080]);
categoryMat=cell2mat(category);
plot([1:11],[categoryMat.NaN]./length(idxTime{1}.Time),'o-')
xticklabels({'13','63','113','163','213','263','313','363','413','463','475'});
ylabel('% Reached Max Current');
xlabel('Square Pulse Time (mSec)');

%% Plot timeseries:
fTime=figure('Position',[0 0 1920 1080]);
plot(Sim{3}.Time, Sim{3}.Data.VSource1);
ylim([0 1.5]);
xlabel('Seconds')
ylabel('Source (V)');
title([num2str(length(Explore{3}{1}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{3}.SimInfo.MaxV) 'V | Random Electrode Placement | Timeseries']);

%% Plot COMM Correlations at Edges:
fCOMMCurr=figure('Position',[0 0 1920 1080]);
count = 1;
for i = 1:length(categories.COMM)
    s(i)=scatter(categories.COMM{i},categories.Curr{i},[],clrs{1}(i,:));
    %     h=lsline; %Linear Fit
    
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
    xlim([1 8])
    ylabel('Current (A)');
    for j = 1:5 %check the first five Explores just for naming purposes
        if ~isempty(Explore{i}{j})
            title([num2str(length(Explore{i}{j}.GraphView.currents)) 'nw | ' num2str(length(Sim)) ' Simulations | ' num2str(Sim{j}.SimInfo.MaxV) 'V | Random Electrode Placement']);
        end
    end
    hold on
    % legend([p2,p3, h],'2nd order Polynomial Fit','3rd order Polynomial Fit','Linear Fit');
    [rCorrelation{i}.COMM,pCorrelation{i}.COMM]=corrcoef(netCOMM{i},netCurrs{i});
%     Legend{i}=strcat([num2str(Explore{i}{3}.IndexTime) ' sec']);
end
% legend(Legend)
legend({'Early','Mid','Late'})

% % log10 Current:
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
gscatter([scatterTemp1 scatterTemp2 scatterTemp3]',cell2mat(categories.COMM)',{scatterCat.Early{:},scatterCat.Mid{:},scatterCat.Late{:}}',[clrs2{1}],[],15,'off')
ylabel('Communicability');
xticks([1 2 3 4]);
xticklabels({'Early','Mid','Late'});

%% Watson Strogatz Analysis:
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

%Split colours into bins
rdif=[lightblue(1)-blue(1)];
r2=lightblue(1):(-rdif)/(length(categories.Clust)-1):blue(1);
gdif=[lightblue(2)-blue(2)];
g2=lightblue(2):(-gdif)/(length(categories.Clust)-1):blue(2);
bdif=[lightblue(3)-blue(3)];
b2=lightblue(3):(-bdif)/(length(categories.Clust)-1):blue(3);

clrs2={[r2; g2; b2]'};
for i =1:length(categories.Clust)
    s(i+1)=scatter(mean(categories.Clust{i}),nanmean(categories.PathLength{i}),[],clrs2{1}(i,:));
    hold on
    e=errorbar(mean(categories.Clust{i}),nanmean(categories.PathLength{i}),std(categories.Clust{i}),'horizontal');
    e1=errorbar(mean(categories.Clust{i}),nanmean(categories.PathLength{i}),nanstd(categories.PathLength{i}));
    e.Color=clrs2{1}(i,:);
    e1.Color=clrs2{1}(i,:);
end

l=legend(s,'Watss-Strogatz Distribution (Ordered to Random)','Early Subgraph','Mid Subgraph','Late Subgraph','location','NorthWest');
ylabel('Avg Path Length');
xlabel('Clustering Coefficient');

%% ANOVAS
% COMM
%Make all the same length
FigList = allchild(groot);

maxlength = max(cellfun(@numel, categories.COMM));
COMMtemp = cellfun(@(v) [v, nan(1, maxlength-numel(v))], categories.COMM, 'UniformOutput', false);

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


