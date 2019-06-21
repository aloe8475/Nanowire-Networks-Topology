function [random, ordered]=createRandom_Ordered_Graphs(e500)
savePath='D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis\Network Explore Functions\';

%% Random Graph Analysis:
%Create random graph with same number of vertices as 500nw network, and the
%avg degree of 500nw network.

A=createRandRegGraph(height(e500.Explore.GraphView.Nodes),round(full(mean(e500.Explore.GraphTheory.DEG))));
random.Adj=A;
ranGraph=graph(A);
random.Graph=ranGraph;
%Random Cluster
[random.Ci,random.Q] = community_louvain(A,1);
[random.GlobalClust,random.AvgLocalClust, random.Clust] = clustCoeff(A);
%Random Path Length
random.Path = path_length(A);
random.AvgPath=mean(random.Path);
%Circuit Rank
random.CircuitRank=numedges(ranGraph) - (numnodes(ranGraph) - 1);
% 
% 

%% Ordered Graph Analysis

numConnections=2;
numNodes=500;

G=WattsStrogatz(numNodes,numConnections,0); %create a ordered graph that is connected to its two nearest neighbours

ordered.Adj=adjacency(G);
%Random Cluster
[ordered.Ci,ordered.Q] = community_louvain(ordered.Adj,1);
[ordered.GlobalClust,ordered.AvgLocalClust, ordered.Clust] = clustCoeff(ordered.Adj);
%Random Path Length
ordered.Path = path_length(ordered.Adj);
ordered.AvgPath=mean(ordered.Path);
ordered.Graph=G;
%Circuit Rank
ordered.CircuitRank=numedges(G) - (numnodes(G) - 1);

rndNetworkExist=1;
save([savePath 'Ordered_Random_Graphs_500nw.mat'],'ordered','random','rndNetworkExist');
end 