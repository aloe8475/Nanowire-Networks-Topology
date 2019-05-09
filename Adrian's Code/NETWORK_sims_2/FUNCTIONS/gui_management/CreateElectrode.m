function ElInfo=CreateElectrode(varargin)
if length(varargin)>1
    %electrode args( name vmax vmin noc freq open duty sigtype tottime step
    %... startoff setfreq vset readfreq vread vstartoff readduty setduty)
    ElInfo.Name=varargin{1};
    ElInfo.Vmax=varargin{2};
    ElInfo.Vmin=varargin{3};
    ElInfo.NoC=varargin{4};
    ElInfo.Frequency=varargin{5};
    ElInfo.OpenCheck=varargin{6};
    ElInfo.Duty=varargin{7};
    ElInfo.SigType=varargin{8};
    ElInfo.Time=varargin{9};
    ElInfo.Step=varargin{10};
    ElInfo.TStartOff=varargin{11};
    ElInfo.SetFreq=varargin{12};
    ElInfo.VSet=varargin{13};
    ElInfo.ReadFreq=varargin{14};
    ElInfo.VRead=varargin{15};
    ElInfo.VStartOff=varargin{16};
    ElInfo.ReadDuty=varargin{17};
    ElInfo.SetDuty=varargin{18};
    ElInfo.SType=varargin{19};
    ElInfo.DType=varargin{20};
else
    handles=varargin{1};
    if isequal(handles.SourceType.Value,1)
        ElInfo.SType=1;ElInfo.DType=0;
        DSel=find(contains(handles.ElectrodeList.String,'Source'),1,'last');
        if ~isempty(DSel)
            LasEl=handles.ElectrodeList.String{DSel};
            LastN=sscanf(LasEl,'Source%d');
        else
            LastN=0;
        end
        
        ElInfo.Name=strcat('Source',num2str(LastN+1));
    else
        ElInfo.SType=0;ElInfo.DType=1;
        SSel=find(contains(handles.ElectrodeList.String,'Drain'),1,'last');
        if ~isempty(SSel)
            LasEl=handles.ElectrodeList.String{SSel};
            LastN=sscanf(LasEl,'Drain%d');
        else
            LastN=0;
        end
        ElInfo.Name=strcat('Drain',num2str(LastN+1));
    end
    
    ElInfo.Vmax=str2double(handles.SourceEdit.String);
    ElInfo.Vmin=str2double(handles.SourceMinEdit.String);
    ElInfo.NoC=str2double(handles.NoCEdit.String);
    ElInfo.Frequency=str2double(handles.FrequencyEdit.String);
    ElInfo.OpenCheck=handles.OpenCheck.Value;
    ElInfo.Duty=str2double(handles.DutyEdit.String);
    ElInfo.SigType=handles.SignalPop.String{handles.SignalPop.Value};
    ElInfo.Time=str2double(handles.TimeEdit.String);
    ElInfo.Step=str2double(handles.StepEdit.String);
    ElInfo.TStartOff=str2double(handles.TStartOffEdit.String);
    ElInfo.SetFreq=str2double(handles.SetFreqEdit.String);
    ElInfo.VSet=str2double(handles.VSetEdit.String);
    ElInfo.ReadFreq=str2double(handles.ReadFreqEdit.String);
    ElInfo.VRead=str2double(handles.VReadEdit.String);
    ElInfo.VStartOff=str2double(handles.VStartOffEdit.String);
    ElInfo.ReadDuty=str2double(handles.ReadDutyEdit.String);
    ElInfo.SetDuty=str2double(handles.SetDutyEdit.String);
    
end



