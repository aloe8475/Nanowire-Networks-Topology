function varargout = callSIMULATIONSETTINGS(varargin)
% CALLSIMULATIONSETTINGS MATLAB code for callSIMULATIONSETTINGS.fig
%      CALLSIMULATIONSETTINGS, by itself, creates a new CALLSIMULATIONSETTINGS or raises the existing
%      singleton*.
%
%      H = CALLSIMULATIONSETTINGS returns the handle to a new CALLSIMULATIONSETTINGS or the handle to
%      the existing singleton*.
%
%      CALLSIMULATIONSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALLSIMULATIONSETTINGS.M with the given input arguments.
%
%      CALLSIMULATIONSETTINGS('Property','Value',...) creates a new CALLSIMULATIONSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before callSIMULATIONSETTINGS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to callSIMULATIONSETTINGS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help callSIMULATIONSETTINGS

% Last Modified by GUIDE v2.5 07-Nov-2018 11:22:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @callSIMULATIONSETTINGS_OpeningFcn, ...
    'gui_OutputFcn',  @callSIMULATIONSETTINGS_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before callSIMULATIONSETTINGS is made visible.
function callSIMULATIONSETTINGS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to callSIMULATIONSETTINGS (see VARARGIN)

% Choose default command line output for callSIMULATIONSETTINGS
handles.output = hObject;

handles=fillSIMULATIONSETTINGS(0,varargin{1},handles);
%electrode name vmax vmin noc freq open duty sigtype tottime step
if isempty(handles.ElectrodesInfo)
    handles.ElectrodesInfo{1}=CreateElectrode('Source1',1,1e-3,1,0,0,25,'Constant',1,0.01...
        ,0,10,5,1,5,0,5,5,1,0);
    handles.ElectrodesInfo{2}=CreateElectrode('Drain1',0,1e-3,1,0,0,25,'Constant',1,0.01...
        ,0,10,5,1,5,0,5,5,0,1);
    handles.ElectrodeList.String={'Source1';'Drain1'};
    handles.ElectrodeList.Value=1;
else
    handles.ElectrodeList.String=cellfun(@(c) c.Name,handles.ElectrodesInfo,'UniformOutput',false);
    handles.ElectrodeList.Value=1;
end

% Update handles structure
guidata(hObject, handles);

set(handles.figure1,'WindowStyle','Modal');

% UIWAIT makes callSIMULATIONSETTINGS wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = callSIMULATIONSETTINGS_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isvalid(hObject)
    varargout{1} = []; %case you click close on the modal window
else
    varargout{1}=handles.SimulationSettings;
    delete(handles.figure1);
end


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimulationSettings=fillSIMULATIONSETTINGS(1,[],handles);

guidata(hObject,handles);

uiresume(handles.figure1);


% --- Executes on button press in PreviewButton.
function PreviewButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for i=1:length(handles.ElectrodesInfo)
[t,v,~]=SignalGenerator(handles.ElectrodesInfo{i});
plot(handles.axes1,t,v);
hold on;
end
XLim=handles.axes1.XLim;
YLim=handles.axes1.YLim;
handles.axes1.XLim=[XLim(1)-XLim(1)*0.05 XLim(2)+XLim(2)*0.05];
handles.axes1.YLim=[YLim(1)-YLim(1)*0.05 YLim(2)+YLim(2)*0.05];
hold off;
xlabel('t(s)');
ylabel('Vsource (v)');
guidata(hObject,handles);


% --- Executes on selection change in ElectrodeList.
function ElectrodeList_Callback(hObject, eventdata, handles)
% hObject    handle to ElectrodeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ElectrodeList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ElectrodeList
Electrode=handles.ElectrodesInfo{handles.ElectrodeList.Value};
handles.SourceEdit.String=num2str(Electrode.Vmax);
handles.SourceMinEdit.String=num2str(Electrode.Vmin);
handles.NoCEdit.String=num2str(Electrode.NoC);
handles.FrequencyEdit.String=num2str(Electrode.Frequency);
handles.OpenCheck.Value=Electrode.OpenCheck;
handles.DutyEdit.String=Electrode.Duty;
handles.SignalPop.Value=find(strcmp(handles.SignalPop.String,Electrode.SigType));
handles.NameText.String=Electrode.Name;
handles.TotalText.String=num2str(Electrode.Time);
handles.TStartOffEdit.String=num2str(Electrode.TStartOff);
handles.SourceType.Value=Electrode.SType;
handles.DrainType.Value=Electrode.DType;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ElectrodeList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ElectrodeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddButton.
function AddButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NewElectrode=CreateElectrode(handles);
[handles.ElectrodesInfo,handles.ElectrodeList.String,...
    handles.ElectrodeList.Value]=...
    AddDeleteModifyElectrodeList('add',handles.ElectrodesInfo,NewElectrode,0);
guidata(hObject,handles);
% --- Executes on button press in RemoveButton.
function RemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelElectrode=handles.ElectrodeList.String{handles.ElectrodeList.Value};
[handles.ElectrodesInfo,handles.ElectrodeList.String,...
    handles.ElectrodeList.Value]=...
    AddDeleteModifyElectrodeList('delete',handles.ElectrodesInfo,SelElectrode,...
    handles.ElectrodeList.Value);
guidata(hObject,handles);


% --- Executes on button press in OpenCheck.
function OpenCheck_Callback(hObject, eventdata, handles)
% hObject    handle to OpenCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OpenCheck


% --- Executes on button press in ModifyButton.
function ModifyButton_Callback(hObject, eventdata, handles)
% hObject    handle to ModifyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelName=handles.ElectrodeList.String{handles.ElectrodeList.Value};
SelElectrode=CreateElectrode(handles);
SelElectrode.Name=SelName;
[handles.ElectrodesInfo,handles.ElectrodeList.String,...
    handles.ElectrodeList.Value]=...
    AddDeleteModifyElectrodeList('modify',handles.ElectrodesInfo,SelElectrode,...
    handles.ElectrodeList.Value);
guidata(hObject,handles);


% --- Executes on selection change in SignalPop.
function SignalPop_Callback(hObject, eventdata, handles)
% hObject    handle to SignalPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SignalPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SignalPop


% --- Executes on button press in UpAllButton.
function UpAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
El=handles.ElectrodesInfo;
for i=1:length(El)
    El{i}.Time=str2double(handles.TimeEdit.String);
end
handles.ElectrodesInfo=El;
guidata(hObject,handles);
