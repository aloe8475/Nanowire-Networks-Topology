function SelSims=Convert_Zdenka_to_Adrian(SelSims,snapshots,SimulationOptions,Connectivity,Components,Stimulus)
%ASK RUOMIN TO HELP CHANGE FOR MULTI ELECTRODES

%Junction Values
    SelSims.Time= Stimulus{1}.TimeAxis;
    SelSims.Data.JunctionRmat=SelSims.Data.JunctionResistance;
    SelSims.Data.JunctionVoltage=SelSims.Data.JunctionVoltages;

%Wire Values
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
    SelSims.Data.VSource1=SelSims.Data.WireVoltages(SelSims.Electrodes(1).PosIndex);
    SelSims.Data.VDrain1=SelSims.Data.WireVoltages(SelSims.Electrodes(2).PosIndex);
    SelSims.Data.ISource1=SelSims.Data.ElectrodeCurrents(:,1); 
    SelSims.Data.IDrain1=SelSims.Data.ElectrodeCurrents(:,2); 

%Adj Matrix & Wire Positions
SelSims.SelLayout.AdjMat=Connectivity.weights;
connLength=length(Connectivity.WireEnds);
connLengthNodes=length(Connectivity.VertexPosition);
xa=Connectivity.WireEnds(:,1);
ya=Connectivity.WireEnds(:,2);
xb=Connectivity.WireEnds(:,3);
yb=Connectivity.WireEnds(:,4);
xa=diag(xa);
xb=diag(xb);
ya=diag(ya);
yb=diag(yb);
xc=(xa+xb)./2; % NEED HELP WITH THIS
yc=(ya+yb)./2; % NEED HELP WITH THIS
% xc=Connectivity.VertexPosition(:,1);
% yc=Connectivity.VertexPosition(:,2);
xi=Connectivity.EdgePosition(1,:);
yi=Connectivity.EdgePosition(2,:);
SelSims.SelLayout.X1=sparse(xa);
SelSims.SelLayout.X2=sparse(xb);
SelSims.SelLayout.Y1=sparse(ya);
SelSims.SelLayout.Y2=sparse(yb);
SelSims.SelLayout.CY=sparse(yc);
SelSims.SelLayout.CX=sparse(xc);

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
SelSims.Settings.Vmax = Stimulus{1}.AmplitudeOn;
SelSims.Settings.Vmin = Stimulus{1}.AmplitudeOff;
if strcmp(Stimulus{1}.BiasType,'AlonPulse')
SelSims.Settings.NoC = Stimulus{1}.NumPulse1+Stimulus.NumPulse2;
SelSims.Settings.SetFreq = Stimulus{1}.Period;

end 
SelSims.Settings.Model = 'Zdenka';
SelSims.Settings.Name  = Connectivity.filename;
SelSims.SimInfo.MaxI=max(max([SelSims.Data.Currents{:}]));
SelSims.SimInfo.MinI=min(min([SelSims.Data.Currents{:}]));
SelSims.SimInfo.MaxV=max(max(SelSims.Data.JunctionVoltage));
SelSims.SimInfo.MinV=min(min(SelSims.Data.JunctionVoltage));
end
