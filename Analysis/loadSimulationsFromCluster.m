%Loading 100 different simulations from cluster

function [Explore,threshold,analysis_type,Sim]=loadSimulationsFromCluster()

loadpath = '/import/silo2/aloe8475/Documents/CODE/Data/Explore Analysis/Time Delay Analysis/';
for i=1:100
    temp=load([loadpath 'ExploreAnalysis_TimeDelay_' num2str(i) '.mat']);
    a=temp.Explore(i,:);
    b=temp.threshold{i};
    for j = 1:length(a)
    Explore{i}{j}=a{j}{i}; %the Explore variable was saved weirdly, so we do this to fix it
    threshold{i}{j}=b{j}{i};
    end 
    Sim{i}=temp.Sim{i};
    analysis_type{i}=temp.analysis_type{i};
    clear temp a b c
end
end

