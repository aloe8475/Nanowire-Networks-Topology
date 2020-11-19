% numInterModuleConnections=[1,2,5,10,15,20,25,50,75,100];
alpha=linspace(0.01,1,20);
% for j = 1:60
for j =1:10
    randSeed(j)=randi(1000);
end
progressbar(0,0)
for i = 1:10
%     count=numInterModuleConnections(i);
    progressbar([],0)
    for j =1:length(alpha)
        %HMN(1,nodes in each module, number of modules, density (lower =
        %sparse),number of levels, probability of rewiring,seed)
        network1{i,j}=HMN(1,100,3,alpha(j),1,0.9,randSeed(i));
        network2{i,j}=HMN(1,60,5,alpha(j),1,0.9,randSeed(i));
        network3{i,j}=HMN(1,30,10,alpha(j),1,0.9,randSeed(i));%,HMN(2,4,3,4,4,0,p(j)),HMN(2,1,7,4,3,0,p(j)),HMN(2,5,4,4,3,0,p(j)),HMN(2,5,2,4,6,0,p(j)),HMN(2,1,7,4,3,0,p(j)),HMN(2,4,9,4,2,0,p(j)),HMN(2,10,2,4,5,0,p(j)),HMN(2,9,2,4,5,0,p(j)),HMN(2,20,2,4,4,0,p(j))};
        progressbar([],j/10)
    end
    progressbar(i/10)
end

networks=[network1;network2;network3];
    % end 
% output=HMN(2,4,3,4,4,0,1);
count=1;
for i = 1:size(networks,1)
    for j = 1:size(networks,2)
%     for j = 1:length(networks{i})
            g=graph(networks{i,j});
            adj=full(adjacency(g));
            mat2np(adj,strjoin(['C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\HierarchicalGraphs\hierarchical_' string(count) '_weighted.pkl'],''),'int8');
            count=count+1;
    end
end 
% clear output
% save(strjoin(['C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\HierarchicalGraphs\hierarchical_' string(length(adj)) '_Nodes.mat'],''))