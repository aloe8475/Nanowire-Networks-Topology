function SelSims=Convert_Zdenka_to_Adrian(SelSims,snapshots,SimulationOptions,Connectivity,Components,Stimulus)

%Junction Values
for i = 1:length(snapshots)
    SelSims.Time= snapshots{i}.Timestamp;
    SelSims.Data.JunctionRmat{i}=snapshots{i}.Resistance;
    SelSims.Data.JunctionCurrents{i}=snapshots{i}.Current;
    SelSims.Data.JunctionVoltage{i}=snapshots{i}.Voltage;
end

%Wire Values
for j = 1:length(SelSims.Data.WireVoltages) %for each timestep
    countSource=0;
    countDrain=0;
    for k = 1:length(SimulationOptions.ContactNodes)
        SelSims.Electrodes(k).PosIndex=SimulationOptions.ContactNodes(k); %this assumes the nodes are the same for each simulation
        if k ==1 %ASK RUOMIN TO CHANGE FOR MULTI ELECTRODES
            SelSims.Electrodes(k).Name=['Source' num2str(countSource)];
        elseif k==2
            SelSims.Electrodes(k).Name=['Drain' num2str(countDrain)];
        end
    end
    SelSims.Data.VSource1(j)=SelSims.Data.WireVoltages{j}(SelSims.Electrodes(1).PosIndex);
    SelSims.Data.VDrain1(j)=SelSims.Data.WireVoltages{j}(SelSims.Electrodes(2).PosIndex);
    SelSims.Data.ISource1(j)=SelSims.Data.WireCurrents{j}{SelSims.Electrodes(1).PosIndex};
    SelSims.Data.IDrain1(j)=SelSims.Data.WireCurrents{j}{SelSims.Electrodes(2).PosIndex};
end

%Adj Matrix & Wire Positions
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

%Settings:
SelSims.Settings.IniRes = Components.offResistance;
SelSims.Settings.Time = SimulationOptions.T;
SelSims.Settings.Step = SimulationOptions.dt;
SelSims.Settings. Ron = Components.onResistance;
SelSims.Settings.Roff = Components.offResistance;
SelSims.Settings.MaxW = Components.maxFlux;
%Ask Joel/Ruomin
%                SelSims.Settings.Tau = 1
%             SelSims.Settings.SigmaW = 0.1000
%         SelSims.Settings.SigmaNoise = 0.1000
%              SelSims.Settings.Alpha = 1.0000e-03
SelSims.Settings.IniW= Components.criticalFlux;
SelSims.Settings.Vmax = Stimulus.AmplitudeOn;
SelSims.Settings.Vmin = Stimulus.AmplitudeOff;
SelSims.Settings.NoC = Stimulus.NumPulse1+Stimulus.NumPulse2;
SelSims.Settings.SetFreq = Stimulus.Period;
SelSims.Settings.Model = 'Zdenka';
SelSims.Settings.Name  = Connectivity.filename;
end
