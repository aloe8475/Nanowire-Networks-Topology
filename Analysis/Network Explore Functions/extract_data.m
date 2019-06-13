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
    if ismember(['IDrain' num2str(i)],fieldnames(network.Simulations{simNum}.Data)) %if we have drains
        save_drain='y';
        IDrain(:,i)=network.Simulations{simNum}.Data.(['IDrain' num2str(i)]);
    end
end

computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        path='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Python Data\';
        % path='D:\alon_\Research\POSTGRAD\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only\Python Data\';
end

SelSims=network.Simulations{simNum};
AdjMat=SelSims.SelLayout.AdjMat;
Time=SelSims.Settings.Time;
if save_drain=='y'
    save_name=[network.Name SelSims.Settings.Model '_' SelSims.Settings.SigType '_' num2str(length(SelSims)) 'SimsOnly_' num2str(SelSims.Settings.Time) '_Sec_' num2str(height(SelSims.Electrodes)) 'Electrodes_Vmax_' num2str(SelSims.Settings.Vmax) '_' network.Simulations{simNum}.Type '_' num2str(simNum) '_' date '_ForPython'];
    save_name=strrep(save_name,':','_');
    save_name=strrep(save_name,'/','_');
    save([path save_name '.mat'],'IDrain','VSource','AdjMat','Time'); %only save data if we have at least 1 source and 1 drain
end
end