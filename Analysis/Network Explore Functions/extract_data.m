function extract_data(network,simNum)
%Training Data:
IDrain1=network.Simulations{simNum}.Data.IDrain1;
IDrain2=network.Simulations{simNum}.Data.IDrain2;
VSource1=network.Simulations{simNum}.Data.VSource1;
VSource2=network.Simulations{simNum}.Data.VSource2;

path='D:\alon_\Research\POSTGRAD\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Python Data\';

save([path network.Simulations{simNum}.Type '_' num2str(simNum) '_ForPython.mat'],'IDrain1','IDrain2','VSource1','VSource2');
end 