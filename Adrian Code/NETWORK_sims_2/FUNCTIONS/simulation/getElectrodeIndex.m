function [elindex,elvalues]=getElectrodeIndex(electrodes)


ElectrodeList=electrodes;

elindex=zeros(1,length(ElectrodeList));
elvalues=elindex;
for i=1:length(elindex)
    elindex(i)=ElectrodeList{i}.IndexNanowire(1);
    elvalues(i)=ElectrodeList{i}.Value;
end


end