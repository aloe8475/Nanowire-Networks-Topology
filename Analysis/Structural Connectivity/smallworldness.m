load('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN.mat')
load('C:\Users\61424\Documents\GitHub\CODE\Data\Organic Networks Connectomes\celegans277neurons.mat')
load('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300_WS.mat')

for i=1:length(AdjMat)
smallworld(i)=small_world_propensity(double(AdjMat{i}));
smallworld_grid(i)=small_world_propensity(double(AdjMat_Grid{i}));
smallworld_random(i)=small_world_propensity(double(AdjMat_Random{i}));
end
cElegansSW=small_world_propensity(celegans277matrix);

save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300nwASN_smallworld.mat','smallworld')
save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\cElegans_smallworld.mat','cElegansSW')
save('C:\Users\61424\Documents\GitHub\CODE\Analysis\Functional Connectivity\Functional Tasks\300_WS_smallworld.mat','smallworld_grid','smallworld_random')


