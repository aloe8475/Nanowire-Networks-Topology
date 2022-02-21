% numInterModuleConnections=[1,2,5,10,15,20,25,50,75,100];
% alpha=linspace(0.01,1,20);
% for j = 1:60

%Load Preconstructed Adj Matrices:
dataLoc='C:/Users/61424/Documents/GitHub/CODE/Data/Functional Connectivity/HierarchicalGraphs/';
mat = dir([dataLoc '*.mat']); 
for q = 1:length(mat)
    cont = load(mat(q).name);
    temp{q}=cont.adj_matrix;
end 

%randomly select 10 networks: These will be used as the module in level 1
%of the hierarchical networks
idx=randperm(40,10);
adj_mats=temp(idx);

%Initialize Parameters for Hierarchical Networks
for j =1:10
    randSeed(j)=randi(1000);
end
numNodes=75;
numModules=2;
numLevels=3;
nModConnections=25;
nLevelConnections=8;
for i = 1:10
%     count=numInterModuleConnections(i);
%     progressbar([],0)
%     for j =1:length(alpha)

    %HMN(1,nodes in each module, number of modules, density (lower =
    %sparse),number of levels, probability of rewiring,seed)
    
%     network1{i}=HMN(1,75,2,0.05,2,0.1,randSeed(i));

    %HMN(2,nodes in each module, number of modules, number of connections in each level
    %number of connections between each level, number of levels, probability of rewiring,module adjMat,seed)
    
    network1{i}=HMN(2, numNodes, numModules, nModConnections,nLevelConnections, numLevels, 1,adj_mats{i}, randSeed(i));

%         network2{i,j}=HMN(1,60,5,alpha(j),1,0.9,randSeed(i));
%         progressbar([],j/10)
%     end
%     progressbar(i/10)
end

networks=[network1];%;network2;network3];
    % end 
% output=HMN(2,4,3,4,4,0,1);
count=1;
for i = 1:size(networks,2)
%     for j = 1:size(networks,2)
%     for j = 1:length(networks{i})
    g=graph(networks{i});
    adj=full(adjacency(g));
    mat2np(adj,strjoin(['C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\HierarchicalGraphs\hierarchical_' string(count) '_' string(numModules) 'x' string(numLevels) '_' string(nModConnections) 'Connections_' string(numNodes) 'Nodes.pkl'],''),'int8');
    count=count+1;
%     end
end 
% clear output
% save(strjoin(['C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\HierarchicalGraphs\hierarchical_' string(length(adj)) '_Nodes.mat'],''))