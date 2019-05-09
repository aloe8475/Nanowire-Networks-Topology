function str=DisplayInfoText(Sim)
TotalTime=max(Sim.Time);
NoNod=nnz(Sim.SelLayout.AdjMat);
NoWires=length(Sim.SelLayout.AdjMat);
NoElectrodes=height(Sim.Electrodes);
MaxV=Sim.SimInfo.MaxV;MinV=Sim.SimInfo.MinV;
MaxI=Sim.SimInfo.MaxI;MinI=Sim.SimInfo.MinI;

str=cell(5,1);
str{1}=strcat('Name:' ,Sim.Name);
str{2}=strcat('NoNodes:',num2str(NoNod),' NoWires:',num2str(NoWires),' NoElectrodes',num2str(NoElectrodes));
str{3}=strcat('TimeIndex:','1'); 
str{4}=strcat('Vmax:',num2str(MaxV,'%.4g'),' Vmin:',num2str(MinV,'%.4g'));
str{5}=strcat('Imax:',num2str(MaxI,'%.4g'),' Imin:',num2str(MinI,'%.4g'));
str{6}=strcat('TotalSimTime:',num2str(TotalTime));

end