
% computer=getenv('computername');
% switch computer
%     case 'W4PT80T2' %if on desktop at uni - Alon
%         savePath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Network Explore Functions\';
%     case '' %if on linux
%         savePath='/suphys/aloe8475/Documents/CODE/Analysis/NetworkExploreFunctions/';
%     case 'LAPTOP-S1BV3HR7'
%         savePath='D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis\Network Explore Functions\';
%         %case '' %--- Add other computer paths (e.g. Mike)
% end

%% Random Graph Analysis:
%Create random graph with same number of vertices as network, and the
%avg degree of network.
% dbstop if error

addpath(genpath('../../'));

        fprintf('Loading Files... \n'); 
load('AgNWGraphTheory.mat','AgNW');

% for i = 1:100
%     temp=round(full(mean(network.Explore.GraphTheory.DEG)));
%     if mod(temp,2)~=0
%         temp=temp+1;
%     end
% A=createRandRegGraph(height(network.Explore.GraphView.Nodes),temp);%this only works with an even number of edges*nodes, so we need to make it even above

%Make 50 versions of WattsStrogatz graphs

if exist('WattsStrogatz_Bootstrap.mat','file')
    load('WattsStrogatz_Bootstrap.mat');
    fprintf('Files Loaded \n'); 
else
    fprintf('Bootstrapping... \n');
    for graph=1:50 %bootstrapping
        %     fprintf([num2str(graph) '\n']);
        Parameters(graph).seed=RandStream.create('mrg32k3a','Seed',randi([1 100])); % set random seed generator
        RandStream.setGlobalStream(Parameters(graph).seed);
        count=0;
        for b = 0:0.1:1
            count = count+1;
            ws{graph}{1,count} = WattsStrogatz(100,round(mean(full(AgNW.Degree{1})/2)),b);
            ws{graph}{2,count} = WattsStrogatz(498,round(mean(full(AgNW.Degree{2})/2)),b);
            ws{graph}{3,count} = WattsStrogatz(1000,round(mean(full(AgNW.Degree{3})/2)),b);
            ws{graph}{4,count} = WattsStrogatz(2000,round(mean(full(AgNW.Degree{4})/2)),b);
        end
        save('WattsStrogatz_Bootstrap.mat','ws','Parameters');
    end
            fprintf('Bootstrapping Completed \n');
end


randomOrdered=cell(50);

% if ~exist('pool','var')
pool=parpool(10);
% end 
fprintf('Starting Graph Theory Analysis\n');

parfor graph=1:50
    for i =1:size(ws{graph},1)
        fprintf([num2str(graph) '\n'])
        for j = 1:size(ws{graph},2)
            A=adjacency(ws{graph}{i,j});
            randomOrdered{graph}.Adj{i,j}=A;
            ranGraph=ws{graph}{i,j};
            randomOrdered{graph}.Graph{i,j}=ranGraph;
            
            %Random Cluster
            [randomOrdered{graph}.Ci{i,j},randomOrdered{graph}.Q{i,j}] = community_louvain(A,1);
            [randomOrdered{graph}.GlobalClust{i,j},randomOrdered{graph}.AvgLocalClust{i,j}, randomOrdered{graph}.Clust] = clustCoeff(A);
            fprintf('Clustering Completed \n');
            
            %Participation Coefficient & Module z-Score
            randomOrdered{graph}.P{i,j} = participation_coef(A,randomOrdered{graph}.Ci{i,j});
            randomOrdered{graph}.MZ{i,j} = module_degree_zscore(A, randomOrdered{graph}.Ci{i,j});
            fprintf('PC & MZ Completed \n');
            
            %Betweenness Centrality:
%             [randomOrdered{graph}.BC{i,j}, randomOrdered{graph}.normBC{i,j}]=betweenness_bin(A);
            
            
            %Communicability:
            randomOrdered{graph}.COMM{i,j} = expm(A);
            fprintf('BC & COMM Completed \n');
            
            %Random Path Length
            randomOrdered{graph}.Path{i,j} = path_length(A);
            randomOrdered{graph}.AvgPath{i,j}=mean(randomOrdered{graph}.Path{i,j});
                    fprintf([num2str(graph) '\n'])
                    fprintf([num2str(i) '\n'])
            
            %Small World Propensity:
            randomOrdered{graph}.SmallWorldProp{i,j}=small_world_propensity(A);
            %Circuit Rank
            randomOrdered{graph}.CircuitRank{i,j}=numedges(ranGraph) - (numnodes(ranGraph) - 1);
                        fprintf('SmallWorld & Circuit Rank Completed \n');

            %             clear ranGraph
        end
    end
    fprintf([num2str(graph) '\n'])
end
save(['Random_Ordered_Graphs_ALL_networks.mat'],'randomOrdered','-v7.3');

clear j i 

for i = 1:50
        for j = 1:size(ws{i},1)
            for k = 1:length(ws{i})
        tempPath(i,j,k)=randomOrdered{i}.AvgPath{j,k};
        tempClust(i,j,k)=randomOrdered{i}.GlobalClust{j,k};
        tempP(i,j,k)=randomOrdered{i}.P{j,k};
        tempMZ(i,j,k)=randomOrdered{i}.MZ{j,k};
        tempBC(i,j,k)=randomOrdered{i}.BC{j,k};
        tempSmallWorld(i,j,k)=randomOrdered{i}.SmallWorldProp{j,k};
        tempCOMM(i,j,k)=randomOrdered{i}.COMM{j,k};
            end 
        end 
end 
        randomOrdered100.AvgPath=squeeze(mean(tempPath));
        randomOrdered100.StdPath(j,:)=squeeze(std(tempPath));
        
       randomOrdered100.AvgGlobalCluster=squeeze(mean(tempClust));
        randomOrdered100.StdGlobalCluster(j,:)=squeeze(std(tempClust));
        
        randomOrdered100.AvgPCoeff(j,:)=squeeze(mean(tempP));
        randomOrdered100.StdPCoeff(j,:)=squeeze(std(tempP));
        
     randomOrdered100.AvgMZ(j,:)=squeeze(mean(tempMZ));
        randomOrdered100.StdMZ(j,:)=squeeze(std(tempMZ));
        
        randomOrdered100.AvgBC(j,:)=squeeze(mean(tempBC));
        randomOrdered100.StdBC(j,:)=squeeze(std(tempBC));
        
        randomOrdered100.AvgSmallWorldProp(j,:)=squeeze(mean(tempSmallWorld));
        randomOrdered100.StdSmallWorldProp(j,:)=squeeze(std(tempSmallWorld));
        
        randomOrdered100.AvgCOMM(j,:)=squeeze(mean(tempCOMM));
        randomOrdered100.StdCOMM(j,:)=squeeze(std(tempCOMM));

save(['Random_Ordered_Graphs_ALL_networks.mat'],'randomOrdered','randomOrdered100','-v7.3');

% %% Ordered Graph Analysis
% for j = 1:100
% numConnections=2;
% numNodes=height(network.Explore.GraphView.Nodes); %same as 500nw network
%
% G=WattsStrogatz(numNodes,numConnections,0); %create a ordered graph that is connected to its two nearest neighbours
%
% ordered{i,j}.Adj=adjacency(G);
% B=ordered{i,j}.Adj;
% %Random Cluster
% [ordered{i,j}.Ci,ordered{i,j}.Q] = community_louvain(B,1);
% [ordered{i,j}.GlobalClust,ordered{i,j}.AvgLocalClust, ordered{i,j}.Clust] = clustCoeff(B);
%
% %Participation Coefficient & Module z-Score
% ordered{i,j}.P = participation_coef(B,ordered{i,j}.Ci);
% ordered{i,j}.MZ = module_degree_zscore(B, ordered{i,j}.Ci);
%
% %Communicability:
% ordered{i,j}.COMM = expm(B);
%
% %Betweenness Centrality:
% [ordered{i,j}.BC, ordered{i,j}.normBC]=betweenness_bin(B);
%
% %Random Path Length
% ordered{i,j}.Path = path_length(B);
% ordered{i,j}.AvgPath=mean(ordered{i,j}.Path);
% ordered{i,j}.Graph=G;
% %Circuit Rank
% ordered{i,j}.CircuitRank=numedges(G) - (numnodes(G) - 1);
% ordered{i,j}.SmallWorldProp=small_world_propensity(B);
% clear G B
% fprintf([num2str{i,j} '\n'])
% end
%
% ordered100.AvgPath=mean([ordered(:).AvgPath]);
% ordered100.StdPath=std([ordered(:).AvgPath]);
%
% ordered100.AvgGlobalClust=mean([ordered(:).GlobalClust]);
% ordered100.StdGlobalClust=std([ordered(:).GlobalClust]);
%
% ordered100.AvgPCoeff=mean([ordered(:).P]);
% ordered100.StdPCoeff=std([ordered(:).P]);
%
% ordered100.AvgMZ=mean([ordered(:).MZ]);
% ordered100.StdMZ=std([ordered(:).MZ]);
%
% ordered100.AvgBC=mean([ordered(:).BC]);
% ordered100.StdBC=std([ordered(:).BC]);
%
% ordered100.AvgSmallWorldProp=mean([ordered(:).SmallWorldProp]);
% ordered100.StdSmallWorldProp=std([ordered(:).SmallWorldProp]);
%
% ordered100.AvgCOMM=mean([ordered(:).COMM]);
% ordered100.StdCOMM=std([ordered(:).COMM]);
%
% save([savePath 'Ordered_Graphs_' num2str(sizeNetwork) 'nw.mat'],'ordered','ordered100','-v7.3');
% end