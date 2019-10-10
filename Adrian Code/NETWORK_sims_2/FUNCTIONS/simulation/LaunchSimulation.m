function [DataSim,NewName]=LaunchSimulation(SimSets)
%% prepare sim

Electrodes=SimSets.Electrodes;
t=SimSets.Time;

[Sim]=fillRfromIV('ini',SimSets);
IVStruct=KirchoffSolver(table2struct(Sim),Electrodes,1);
Sim=UpdateSimTable(table2struct(Sim),IVStruct,Electrodes,0);


tic;
for i=2:length(t)
    RStruct=fillRfromIV('loop',Sim(i-1,:),SimSets.Settings);
    IVStruct=KirchoffSolver(RStruct,Electrodes,i);
    disp(strcat('SumRule: ',num2str(IVStruct.SumRule)));
    disp(strcat('Sim Iteration:',num2str(i),'_Iterations left :',num2str(length(t)-i)));
    Sim(i,:)=UpdateSimTable(RStruct,IVStruct,Electrodes,t(i));
    [Electrodes,AllFlag]=DetectThresoldCurrent(Sim,Electrodes,SimSets);
    if isequal(1,AllFlag)
        SimSets.Time=t(1:i);
        SimSets.Settings.Time=t(i);
        t=t(1:i);
        break;
    end
        
end
toc;


numel=length(SimSets.Settings.ElectrodesInfo);

ElectrodePosStr =num2str([SimSets.Electrodes(:,:).PosIndex]');
ElectrodePosStr=regexprep(ElectrodePosStr, '  ', '_'); %Alon to fix 07/06/19

NewName=[SimSets.Settings.PreName '_M_' SimSets.Settings.Model '_Type_' SimSets.Settings.SigType ...
    '_' num2str(numel) '_Electrodes_Positions_' ElectrodePosStr '_Vmax_' num2str(SimSets.Settings.Vmax) 'V' ,...
    '_T_' num2str(length(t)) '_Created_' date];

DataSim=Sim;

end