function [randomOrdered, randomOrdered100, Parameters]=createRandom_Ordered_Graphs(AgNW)
computer=getenv('computername');
switch computer
    case 'W4PT80T2' %if on desktop at uni - Alon
        savePath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Network Explore Functions\';
    case '' %if on linux
        savePath='/suphys/aloe8475/Documents/CODE/Analysis/NetworkExploreFunctions/';
    case 'LAPTOP-S1BV3HR7'
        savePath='D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis\Network Explore Functions\';
        %case '' %--- Add other computer paths (e.g. Mike)
end

Parameters.seed=RandStream.create('mrg32k3a','Seed',42); % set random seed generator
RandStream.setGlobalStream(Parameters.seed);

%% Random Graph Analysis:
%Create random graph with same number of vertices as network, and the
%avg degree of network.
dbstop if error




% for i = 1:100
%     temp=round(full(mean(network.Explore.GraphTheory.DEG)));
%     if mod(temp,2)~=0
%         temp=temp+1;
%     end
% A=createRandRegGraph(height(network.Explore.GraphView.Nodes),temp);%this only works with an even number of edges*nodes, so we need to make it even above

%Make 50 versions of WattsStrogatz graphs

for graph=1%:50 %bootstrapping
    count=0;
    for b = 0:0.1:1
        count = count+1;
        ws{graph}{1,count} = WattsStrogatz(100,round(mean(full(AgNW.Degree{1})/2)),b);
        ws{graph}{2,count} = WattsStrogatz(498,round(mean(full(AgNW.Degree{2})/2)),b);
%         ws{graph}{3,count} = WattsStrogatz(1000,round(mean(full(AgNW.Degree{3})/2)),b);
%         ws{graph}{4,count} = WattsStrogatz(2000,round(mean(full(AgNW.Degree{4})/2)),b);
    end
end
pool=parpool(1);
parfor graph=1:50
    fprintf(num2str(graph))
    for i =1:size(ws,1)
        for j = 1:size(ws,2)
            A=adjacency(ws{graph}{i,j});
            randomOrdered(graph).Adj{i,j}=A;
            ranGraph=ws{graph}{i,j};
            randomOrdered(graph).Graph{i,j}=ranGraph;
            %Random Cluster
            [randomOrdered(graph).Ci{i,j},randomOrdered(graph).Q{i,j}] = community_louvain(A,1);
            [randomOrdered(graph).GlobalClust{i,j},randomOrdered(graph).AvgLocalClust{i,j}, randomOrdered(graph).Clust] = clustCoeff(A);
            %Participation Coefficient & Module z-Score
            randomOrdered(graph).P{i,j} = participation_coef(A,randomOrdered(graph).Ci{i,j});
            randomOrdered(graph).MZ{i,j} = module_degree_zscore(A, randomOrdered(graph).Ci{i,j});
            %Betweenness Centrality:
            [randomOrdered(graph).BC{i,j}, randomOrdered(graph).normBC{i,j}]=betweenness_bin(A);
            
            %Communicability:
            randomOrdered(graph).COMM{i,j} = expm(A);
            
            %Random Path Length
            randomOrdered(graph).Path{i,j} = path_length(A);
            randomOrdered(graph).AvgPath{i,j}=mean(randomOrdered(graph).Path{i,j});
            
            %Small World Propensity:
            randomOrdered(graph).SmallWorldProp{i,j}=small_world_propensity(A);
            %Circuit Rank
            randomOrdered(graph).CircuitRank{i,j}=numedges(ranGraph) - (numnodes(ranGraph) - 1);
%             clear ranGraph
        end
    end
    fprintf([num2str(graph) '\n'])
end

for i = 1:50
    for j = 1:length(ws{1})
randomOrdered100.AvgPath(j,:)=mean([randomOrdered(i).AvgPath{j,:}]);
randomOrdered100.StdPath(j,:)=std([randomOrdered(i).AvgPath{j,:}]);

randomOrdered100.AvgGlobalClust(j,:)=mean([randomOrdered(i).GlobalClust{j,:}]);
randomOrdered100.StdGlobalClust(j,:)=std([randomOrdered(i).GlobalClust{j,:}]);

randomOrdered100.AvgPCoeff(j,:)=mean([randomOrdered(i).P{j,:}]);
randomOrdered100.StdPCoeff(j,:)=std([randomOrdered(i).P{j,:}]);

randomOrdered100.AvgMZ(j,:) = mean([randomOrdered(i).MZ{j,:}]);
randomOrdered100.StdMZ(j,:) = std([randomOrdered(i).MZ{j,:}]);

randomOrdered100.AvgBC(j,:)=mean([randomOrdered(i).BC{j,:}]);
randomOrdered100.StdBC(j,:)=std([randomOrdered(i).BC{j,:}]);

randomOrdered100.AvgSmallWorldProp(j,:)=mean([randomOrdered(i).SmallWorldProp{j,:}]);
randomOrdered100.StdSmallWorldProp(j,:)=std([randomOrdered(i).SmallWorldProp{j,:}]);

randomOrdered100.AvgCOMM(j,:)=mean([randomOrdered(i).COMM{j,:}]);
randomOrdered100.StdCOMM(j,:)=std([randomOrdered(i).COMM{j,:}]);
    end 
end 
save([savePath 'Random_Ordered_Graphs_ALL_networks.mat'],'randomOrdered','randomOrdered100','-v7.3');
end

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