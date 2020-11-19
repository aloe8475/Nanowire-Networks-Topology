<<<<<<< HEAD
% load('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_multipleNetworks.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Organic Networks Connectomes\celegans277neurons.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\300nwASN_Evolution_multipleNetworks.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\BetaSweepNetworks.mat')
% load('GeneticEvo_BestNetworks.mat')
load('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\ModularGraphs\300NW_9modules_100sets.mat')

count=1;
for i=1:size(AdjMat,1)
    for j=1:size(AdjMat,2)
%     B=cell2mat(AdjMat(i));
%     if length(B) > 1
        if length(AdjMat{i,j})>1
            smallworld(count)=small_world_propensity(logical(AdjMat{i,j}));
        else
            smallworld(count)=nan;
        end
%     else
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
        count=count+1;
    end 
end 
% end
% cElegansSW=small_world_propensity(celegans277matrix);

% save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_multipleNetworks_smallworld.mat','smallworld')
% save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\cElegans_multipleNetworks_smallworld.mat','cElegansSW')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\300_WS_AllASN_smallworld.mat','smallworld_grid','smallworld_random')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\300nwASN_Evolution_multipleNetworks_smallworld.mat','smallworld')
% save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\BetaSweepNetworks_smallworld.mat','smallworld_sweep')
save('C:\Users\61424\Documents\GitHub\CODE\Data\Functional Connectivity\ModularGraphs\300NW_9modules_100sets_smallworld.mat','smallworld')

=======
% load('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_multipleNetworks.mat')
load('C:\Users\61424\Documents\GitHub\CODE\Data\Organic Networks Connectomes\celegans277neurons.mat')
% load('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300_WS.mat')

for i=1:length(AdjMat)
    for j = 1:length(AdjMat(i,:))
        smallworld(i,j)=small_world_propensity(double(AdjMat{i,j}));
%         smallworld_grid(i)=small_world_propensity(double(AdjMat_Grid{i,j}));
%         smallworld_random(i)=small_world_propensity(double(AdjMat_Random{i}));
    end 
end
cElegansSW=small_world_propensity(celegans277matrix);

% save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_multipleNetworks_smallworld.mat','smallworld')
save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\cElegans_multipleNetworks_smallworld.mat','cElegansSW')
% save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300_WS_smallworld.mat','smallworld_grid','smallworld_random')
% 

>>>>>>> ddce18adac73a8f52aef01ad38c0eaa9b5f56778
