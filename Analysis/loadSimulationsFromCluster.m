%Loading 100 different simulations from cluster

function Sims=loadSimulationsFromCluster()

loadpath = 'C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Data\Raw\Simulations\Zdenka\Variable Time Delay\';
for i=1:100
load([loadpath 'SelSims_TimeDelay_' num2str(i) '.mat']);
Sims{i}=SelSims;

end 
end 



