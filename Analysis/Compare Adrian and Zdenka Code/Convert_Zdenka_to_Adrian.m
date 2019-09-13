function SelSims=Convert_Zdenka_to_Adrian(SelSims,snapshots,SimulationOptions)

for i = 1:length(snapshots)
SelSims.Time= snapshots{i}.Timestamp;
SelSims.Data.JunctionRmat{i}=snapshots{i}.Resistance;
SelSims.Data.JunctionCurrents{i}=snapshots{i}.Current;
SelSims.Data.JunctionVoltage{i}=snapshots{i}.Voltage;
end 

for j = 1:length(SelSims.Data.WireVoltages) %for each timestep
SelSims.Electrodes=SimulationOptions.ContactNodes; %this assumes the nodes are the same for each simulation
SelSims.Data.VSource1(j)=SelSims.Data.WireVoltages{j}(SelSims.Electrodes(1));
SelSims.Data.VDrain1(j)=SelSims.Data.WireVoltages{j}(SelSims.Electrodes(2));
SelSims.Data.ISource1(j)=SelSims.Data.WireCurrents{j}{SelSims.Electrodes(1)};
SelSims.Data.IDrain1(j)=SelSims.Data.WireCurrents{j}{SelSims.Electrodes(2)};
end 

SelSims.SelLayout.AdjMat=Connectivity.weights;
connLength=length(Connectivity.WireEnds);
connLengthNodes=length(Connectivity.VertexPosition);
xa=Connectivity.WireEnds(:,1);
ya=Connectivity.WireEnds(:,2);
xb=Connectivity.WireEnds(:,3);
yb=Connectivity.WireEnds(:,4);
xc=Connectivity.VertexPosition(:,1);
yc=Connectivity.VertexPosition(:,2);
xi=Connectivity.EdgePosition(1,:);
yi=Connectivity.EdgePosition(2,:);
SelSims.SelLayout.X1=sparse(xa);
SelSims.SelLayout.X2=sparse(xb);
SelSims.SelLayout.Y1=sparse(ya);
SelSims.SelLayout.Y2=sparse(yb);
end 
