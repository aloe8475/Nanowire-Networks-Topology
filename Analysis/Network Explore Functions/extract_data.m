function extract_data(network,simNum)
%Training Data:
numElec=height(network.Simulations{simNum}.Electrodes);

% IDrain1=network.Simulations{simNum}.Data.IDrain1;
% VSource1=network.Simulations{simNum}.Data.VSource1;
drainIndex=find(contains(network.Simulations{simNum}.Electrodes.Name,'Drain'));
sourceIndex=find(contains(network.Simulations{simNum}.Electrodes.Name,'Source'));

for i = 1:length(sourceIndex)
VSource(:,i)=network.Simulations{simNum}.Data.(['VSource' num2str(i)]);
end 
for i=1:length(drainIndex)
    IDrain(:,i)=network.Simulations{simNum}.Data.(['IDrain' num2str(i)]);
end 

computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        path='C:\Users\aloe8475\Documents\GitHub\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Python Data';
% path='D:\alon_\Research\POSTGRAD\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Python Data\';
end 
save([path network.Simulations{simNum}.Type '_' num2str(simNum) '_ForPython.mat'],'IDrain','VSource');
end 