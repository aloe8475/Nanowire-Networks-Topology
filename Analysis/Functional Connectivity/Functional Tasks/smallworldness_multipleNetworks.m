% load('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_multipleNetworks.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Organic Networks Connectomes\celegans277neurons.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\300nwASN_Evolution_multipleNetworks.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\Betasweep\WattsStrogatz_25deg_NWNMC_allVolts.mat')
load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\Network Comparison\AllNetworks_AdjMats.mat')

% load('GeneticEvo_BestNetworks.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\ModularGraphs\300NW_9modules_NWN_10sets.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\DensityNWNs\VaryingDensity300NWset2_WEIGHTEDSUBGRAPHS_NWN_allVolts.mat')
%%
% %% if networks are 100x300x300 doubles instead of cells:
a=num2cell(AdjMatSparseCrossbar,[2 3]);
for i=1:length(a)
    b{i}=squeeze(a{i});
end 
AdjMatSparseCrossbar=b;
%% 
% smallworld=small_world_propensity(logical(AdjMatNWN),'bin');
% smallworldWS=small_world_propensity(logical(AdjMatWS),'bin');
% smallworldElegans=small_world_propensity(logical(AdjMatElegans),'bin');
% smallworldCrossbar=small_world_propensity(logical(AdjMatCrossbar),'bin');
% smallworldSparseCrossbar=small_world_propensity(logical(AdjMatSparseCrossbar),'bin');
% smallworldER=small_world_propensity(logical(AdjMatER),'bin');
% smallworldRewired=small_world_propensity(logical(AdjMatRewired),'bin');

%%
AdjMats={AdjMatNWN,AdjMatRewired,AdjMatSparseCrossbar};
for i=1:size(AdjMats,2)
    for j = 1:size(AdjMats{2},2)
        smallworld(i,j)=small_world_propensity(logical(AdjMats{i}{j}),'bin'); %for binary networks uncomment
%         smallworld(i)=nan;
%     end
%     for j = 1:5
%         for k = 1:50
%             smallworld_sweep(i,j,k)=small_world_propensity(double(AdjMat_WS{i,j,k}.graphs));
%         temp=logical(AdjMat_Grid{i}.grid);
%         temp2=logical(AdjMat_Random{i}.random);
%         temp3=logical(AdjMat_BA{i}.grid);   
%         smallworld_grid(i)=small_world_propensity(temp);
%         smallworld_random(i)=small_world_propensity(temp2);
%         smallworld_ba(i)=small_world_propensity(temp3);
%         end 
%         count=count+1;
    end 
end 
%%
% end
% cElegansSW=small_world_propensity(celegans277matrix);
% save('C:\Users\61424\Desktop\modularNWs_MultiTasking_AdjMat_smallworld.mat','smallworld')
save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\Network Comparison\AllNetworks_AdjMats_smallworld.mat')

% save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_multipleNetworks_smallworld.mat','smallworld')
% save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\cElegans_multipleNetworks_smallworld.mat','cElegansSW')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\300_WS_AllASN_smallworld.mat','smallworld_grid','smallworld_random')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\300nwASN_Evolution_multipleNetworks_smallworld.mat','smallworld')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\Betasweep\WattsStrogatz_25deg_NWNMC_allVolts_smallworld.mat','smallworld')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\ModularGraphs\300NW_9modules_NWN_10sets_smallworld.mat','smallworld')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\DensityNWNs\VaryingDensity300NWset2_WEIGHTEDSUBGRAPHS_NWN_allVolts_smallworld.mat','smallworld')
