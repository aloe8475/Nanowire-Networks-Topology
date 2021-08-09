%Loading 100 different simulations from cluster

function [Explore,threshold,analysis_type,Sim]=loadSimulationsFromCluster(type,pathLengthType)

switch type
    case 'DC'
loadpath = '/import/silo2/aloe8475/Documents/CODE/Data/Explore Analysis/SingleDCPulse/';
for i=1:100
    temp=load([loadpath 'ExploreAnalysis_DC_' pathLengthType '_' num2str(i) '.mat']);
    a=temp.Explore{i};
    b=temp.threshold{i};
    Explore{i}=a; 
    threshold{i}=b; 
    Sim{i}=temp.Sim{i};
    analysis_type='c';
    clear temp a b c
end
       
    case 'Pulse'
loadpath = '/import/silo2/aloe8475/Documents/CODE/Data/Explore Analysis/DCPulse/';
for i=1:100
    temp=load([loadpath 'ExploreAnalysis_DCandWait_' pathLengthType '_' num2str(i) '.mat']);
    a=temp.Explore{i};
    b=temp.threshold{i};
    for j = 1:length(a)
    Explore{i}{j}=a{j}{1}; %the Explore variable was saved weirdly, so we do this to fix it
    threshold{i}{j}=b{j}{1};
    end 
    Sim{i}=temp.Sim{i};
    analysis_type='e';
    clear temp a b c
end

    case 'Time Delay'
loadpath = '/import/silo2/aloe8475/Documents/CODE/Data/Explore Analysis/Time Delay Analysis/';
% loadpathSim='/import/silo2/aloe8475/Documents/CODE/Data/Raw/Simulations/Zdenka/Variable Time Delay/';
for i=1:100
    temp=load([loadpath 'ExploreAnalysis_TimeDelay_' pathLengthType '_' num2str(i) '.mat']);
%     temp2=load([loadpathSim 'SelSims_TimeDelay_' pathLengthType '_' num2str(i) '.mat']);
    for j = 1:200
    a=temp.Explore{i};
    b=temp.threshold{i};
    for j = 1:length(a)
    Explore{i}{j}=a{j}; %the Explore variable was saved weirdly, so we do this to fix it
    threshold{i}{j}=b{j};
    end 
    Sim{i}=temp.newSim{i};
    analysis_type='t';
    end 
    clear  a b 
end
    cd(loadpath)
end

