function nrt=UpdateSimTable(Rstruct,IVstruct,Electrodes,curr_t)
nrt=table();

nrt.AdjMat{1}=Rstruct.AdjMat;
nrt.Gmat{1}=Rstruct.Gmat;
nrt.Wmat{1}=Rstruct.Wmat;
nrt.Rmat{1}=Rstruct.Rmat;

nrt.SumRule(1)=IVstruct.SumRule;
nrt.Currents{1}=IVstruct.Currents;
nrt.Voltages{1}=IVstruct.Voltages;
Numel=height(Electrodes);
for i=1:Numel
    nrt.(strcat('I',Electrodes.Name{i}))=IVstruct.(strcat('I',Electrodes.Name{i}));
    nrt.(strcat('V',Electrodes.Name{i}))=IVstruct.(strcat('V',Electrodes.Name{i}));
end
nrt.VoltDif{1}=IVstruct.VoltDif;
nrt.time(1)=curr_t;





end