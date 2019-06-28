function [random, ordered, random100, ordered100]=createRandom_Ordered_Graphs(explore)
savePath='D:\alon_\Research\POSTGRAD\PhD\CODE\Analysis\Network Explore Functions\';

%% Random Graph Analysis:
%Create random graph with same number of vertices as 'explore' network, and the
%avg degree of 'explore' network.
dbstop if error

for i = 1:100
A=createRandRegGraph(height(explore.Explore.GraphView.Nodes),round(full(mean(explore.Explore.GraphTheory.DEG))));
random(i).Adj=A;
ranGraph=graph(A);
random(i).Graph=ranGraph;
%Random Cluster
[random(i).Ci,random(i).Q] = community_louvain(A,1);
[random(i).GlobalClust,random(i).AvgLocalClust, random(i).Clust] = clustCoeff(A);
%Participation Coefficient:
random(i).P = participation_coef(A,random(i).Ci);

%Betweenness Centrality:
[random(i).BC, random(i).normBC]=betweenness_bin(A);

%Communicability:
random(i).COMM = expm(A);

%Random Path Length
random(i).Path = path_length(A);
random(i).AvgPath=mean(random(i).Path);

%Small World Propensity:
random(i).SmallWorldProp=small_world_propensity(A);
%Circuit Rank
random(i).CircuitRank=numedges(ranGraph) - (numnodes(ranGraph) - 1);
clear A ranGraph
fprintf([num2str(i) '\n'])
end 

random100.AvgPath=mean([random(:).AvgPath]);
random100.StdPath=std([random(:).AvgPath]);

random100.AvgGlobalClust=mean([random(:).GlobalClust]);
random100.StdGlobalClust=std([random(:).GlobalClust]);


random100.AvgPCoeff=mean([random(:).P]);
random100.StdPCoeff=std([random(:).P]);

random100.AvgBC=mean([random(:).BC]);
random100.StdBC=std([random(:).BC]);

random100.AvgSmallWorldProp=mean([random(:).SmallWorldProp]);
random100.StdSmallWorldProp=std([random(:).SmallWorldProp]);

random100.AvgCOMM=mean([random(:).COMM]);
random100.StdCOMM=std([random(:).COMM]);

save([savePath 'Random_Graphs_500nw.mat'],'random','random100');


%% Ordered Graph Analysis
for j = 1:100
numConnections=2;
numNodes=500;

G=WattsStrogatz(numNodes,numConnections,0); %create a ordered graph that is connected to its two nearest neighbours

ordered(j).Adj=adjacency(G);
B=ordered(j).Adj;
%Random Cluster
[ordered(j).Ci,ordered(j).Q] = community_louvain(B,1);
[ordered(j).GlobalClust,ordered(j).AvgLocalClust, ordered(j).Clust] = clustCoeff(B);

%Participation Coefficient:
ordered(j).P = participation_coef(B,ordered(j).Ci);

%Communicability:
ordered(i).COMM = expm(B);

%Betweenness Centrality:
[ordered(i).BC, ordered(i).normBC]=betweenness_bin(B);

%Random Path Length
ordered(j).Path = path_length(B);
ordered(j).AvgPath=mean(ordered(j).Path);
ordered(j).Graph=G;
%Circuit Rank
ordered(j).CircuitRank=numedges(G) - (numnodes(G) - 1);
ordered(i).SmallWorldProp=small_world_propensity(B);
clear G B
fprintf([num2str(j) '\n'])
end 

ordered100.AvgPath=mean([ordered(:).AvgPath]);
ordered100.StdPath=std([ordered(:).AvgPath]);

ordered100.AvgGlobalClust=mean([ordered(:).GlobalClust]);
ordered100.StdGlobalClust=std([ordered(:).GlobalClust]);

ordered100.AvgPCoeff=mean([ordered(:).P]);
ordered100.StdPCoeff=std([ordered(:).P]);

ordered100.AvgBC=mean([ordered(:).BC]);
ordered100.StdBC=std([ordered(:).BC]);

ordered100.AvgSmallWorldProp=mean([ordered(:).SmallWorldProp]);
ordered100.StdSmallWorldProp=std([ordered(:).SmallWorldProp]);

ordered100.AvgCOMM=mean([ordered(:).COMM]);
ordered100.StdCOMM=std([ordered(:).COMM]);

save([savePath 'Ordered_Graphs_500nw.mat'],'ordered','ordered100');
end 