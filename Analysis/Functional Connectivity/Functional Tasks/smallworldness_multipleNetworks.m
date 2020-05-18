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

