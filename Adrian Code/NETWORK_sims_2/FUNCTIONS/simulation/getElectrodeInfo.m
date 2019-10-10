function [elout]=getElectrodeInfo(network,sim)
%%  - Checks if source and drain are connected
%%  - Chooses the nanowire with the less number of junctions as the electrode
%%  - retrieves the indexes for the nanowires connected to the electrodes
%%  - retrieves the values for source and drain as arrays Vsource, Vdrain.
elout=struct();
ElectrodeList=sim.Electrodes;
Settings=sim.Settings;
LayOut=network.LayOut;
Domains=network.Domains;
Vsource=Settings.Source;
Vdrain=Settings.Drain;

elindex=zeros(1,length(ElectrodeList));
elvalues=elindex;
for i=1:length(elindex)
   
    elout.values(i)=ElectrodeList{i}.Value;
    elout.
end


end