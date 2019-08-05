function varargout=fillSIMULATIONSETTINGS(ini,inset,handles)
%fill simulation settings
%shared function between main_sims and callNSIMULATIONSETTINGS

switch ini
    case 0 %modal figure initial
        handles.PreNameEdit.String=inset.PreName;
        handles.ElectroType=inset.ElectrodeType;
        if isequal(inset.ElectrodeType,'point')
            i1=1;i2=0;
        else
            i1=0;i2=1;
        end
        handles.OneNWButton.Value=i1;
        handles.SquareNWButton.Value=i2;
        
        handles.IniResEdit.String=num2str(inset.IniRes);
        handles.TimeEdit.String=num2str(inset.Time);
        handles.StepEdit.String=num2str(inset.Step);
        handles.RonEdit.String=num2str(inset.Ron);
        handles.RoffEdit.String=num2str(inset.Roff);
        handles.MaxWEdit.String=num2str(inset.MaxW);
        handles.MobilityEdit.String=num2str(inset.Mobility);
        handles.TauEdit.String=num2str(inset.Tau);
        handles.SigmaWEdit.String=num2str(inset.SigmaW);
        handles.SigmaNoiseEdit.String=num2str(inset.SigmaNoise);
        handles.AlphaEdit.String=num2str(inset.Alpha);
        handles.IniWEdit.String=num2str(inset.IniW);
        handles.MaxIEdit.String=num2str(inset.MaxI);
        handles.PowEdit.String=num2str(inset.Pow);
        
        
        handles.ElectrodesInfo=inset.ElectrodesInfo;
        handles.SourceEdit.String=num2str(inset.Vmax);
        handles.SourceMinEdit.String=num2str(inset.Vmin);
        handles.NoCEdit.String=num2str(inset.NoC);
        handles.DutyEdit.String=num2str(inset.Duty);
        handles.FrequencyEdit.String=num2str(inset.Frequency);
        handles.DutyEdit.String=num2str(inset.Duty);
        handles.SignalPop.Value=find(strcmp(handles.SignalPop.String,inset.SigType));
        handles.OpenCheck.Value=inset.OpenCheck;
        handles.NameText.String=inset.Name;
        handles.TStartOffEdit.String=inset.TStartOff;
        handles.SetFreqEdit.String=inset.SetFreq;
        handles.VSetEdit.String=inset.VSet;
        handles.ReadFreqEdit.String=inset.ReadFreq;
        handles.VReadEdit.String=inset.VRead;
        handles.VStartOffEdit.String=inset.VStartOff;
        handles.ReadDutyEdit.String=inset.ReadDuty;
        handles.SetDutyEdit.String=inset.SetDuty;
        handles.SourceType.Value=inset.SType;
        handles.DrainType.Value=inset.DType;
        
        handles.FactorZEdit.String=num2str(inset.FactorZ);
        handles.WFormZEdit.String=num2str(inset.WForm);
        handles.WDissolveZEdit.String=num2str(inset.WDissolve);
        handles.ModelPop.Value=find(strcmp(handles.ModelPop.String,inset.Model));
        
        varargout{1}=handles;
        
    case 1 %filling settings from modal
        sim.PreName=handles.PreNameEdit.String;
        sim.OneNW=handles.OneNWButton.Value;
        sim.SquareNW=handles.SquareNWButton.Value;
        if isequal(sim.OneNW,1)
            aux='point';
        else
            aux='square';
        end
        sim.ElectrodeType=aux;
        sim.IniRes=str2double(handles.IniResEdit.String);
        sim.Time=str2double(handles.TimeEdit.String);
        sim.Step=str2double(handles.StepEdit.String);
        sim.Ron=str2double(handles.RonEdit.String);
        sim.Roff=str2double(handles.RoffEdit.String);
        sim.MaxW=str2double(handles.MaxWEdit.String);
        sim.Mobility=str2double(handles.MobilityEdit.String);
        sim.Tau=str2double(handles.TauEdit.String);
        sim.SigmaW=str2double(handles.SigmaWEdit.String);
        sim.SigmaNoise=str2double(handles.SigmaNoiseEdit.String);
        sim.Alpha=str2double(handles.AlphaEdit.String);
        sim.IniW=str2double(handles.IniWEdit.String);
        sim.Factor=sim.Mobility*sim.Ron/sim.MaxW;
        sim.MaxI=str2double(handles.MaxIEdit.String);
        sim.Pow=str2double(handles.PowEdit.String);
        
        sim.WForm=str2double(handles.WFormZEdit.String);
        sim.FactorZ=str2double(handles.FactorZEdit.String);
        sim.WDissolve=str2double(handles.WDissolveZEdit.String);
        
        sim.ElectrodesInfo=handles.ElectrodesInfo;
        sim.Vmax=str2double(handles.SourceEdit.String);
        sim.Vmin=str2double(handles.SourceMinEdit.String);
        sim.Frequency=str2double(handles.FrequencyEdit.String);
        sim.NoC=str2double(handles.NoCEdit.String);
        sim.Duty=str2double(handles.DutyEdit.String);
        sim.OpenCheck=handles.OpenCheck.Value;
        sim.Name=handles.NameText.String;
        sim.TStartOff=str2double(handles.TStartOffEdit.String);
        sim.SetFreq=str2double(handles.SetFreqEdit.String);
        sim.VSet=str2double(handles.VSetEdit.String);
        sim.VRead=str2double(handles.VReadEdit.String);
        sim.ReadFreq=str2double(handles.ReadFreqEdit.String);
        sim.VStartOff=str2double(handles.VStartOffEdit.String);
        sim.ReadDuty=str2double(handles.ReadDutyEdit.String);
        sim.SetDuty=str2double(handles.SetDutyEdit.String);
        sim.SType=handles.SourceType.Value;
        sim.DType=handles.DrainType.Value;
        
        sim.Model=handles.ModelPop.String{handles.ModelPop.Value};
        sim.SigType=handles.SignalPop.String{handles.SignalPop.Value};
        varargout{1}=sim;
    case 2 %creation of struct in main program
        sim.PreName='Sim';
        sim.ElectrodeType='point';
        sim.OneNW=1;
        sim.SquareNW=0;
        
        sim.IniRes=5e6;
        sim.Time=4;%1;
        sim.Step=0.01;
        sim.Ron=5e3;
        sim.Roff=5e6;
        sim.MaxW=5e-9;
        sim.Mobility=0.5e-12;
        sim.Tau=1;
        sim.SigmaW=0.1;
        sim.SigmaNoise=0.1;
        sim.IniW=1e-11;
        sim.Alpha=1e-3;
        sim.Factor=sim.Mobility*sim.Ron/sim.MaxW;
        sim.MaxI=1e-4;
        sim.WForm=5e-9;
        sim.FactorZ=1;
        sim.WDissolve=1e-9;
        sim.Model='Zdenka';%'HP';
        sim.Pow=0;
        
        
        sim.ElectrodesInfo={};
        sim.Vmax=1;
        sim.Vmin=1e-3;
        sim.NoC=1;
        sim.Frequency=0;
        sim.Name='Source1';
        sim.OpenCheck=0;
        sim.Duty=50;
        sim.SigType='Constant';
        sim.TStartOff=0;
        sim.SetFreq=10;
        sim.VSet=5;
        sim.ReadFreq=1;
        sim.VRead=5;
        sim.VStartOff=0;
        sim.ReadDuty=5;
        sim.SetDuty=5;
        sim.SType=1;
        sim.DType=0;
        
        
        
        varargout{1}=sim;
    otherwise
        return;
end
end