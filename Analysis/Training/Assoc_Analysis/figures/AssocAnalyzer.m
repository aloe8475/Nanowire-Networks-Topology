function varargout = AssocAnalyzer(varargin)
% ASSOCANALYZER MATLAB code for AssocAnalyzer.fig
%      ASSOCANALYZER, by itself, creates a new ASSOCANALYZER or raises the existing
%      singleton*.
%
%      H = ASSOCANALYZER returns the handle to a new ASSOCANALYZER or the handle to
%      the existing singleton*.
%
%      ASSOCANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSOCANALYZER.M with the given input arguments.
%
%      ASSOCANALYZER('Property','Value',...) creates a new ASSOCANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AssocAnalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AssocAnalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AssocAnalyzer

% Last Modified by GUIDE v2.5 27-Jun-2019 15:42:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AssocAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @AssocAnalyzer_OutputFcn, ...
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


% --- Executes just before AssocAnalyzer is made visible.
function AssocAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AssocAnalyzer (see VARARGIN)

% Choose default command line output for AssocAnalyzer
handles.output = hObject;
handles.refList={};
handles.defColormap=jet;
handles.Scaling='direct';
handles.PlotType='PAll';
handles.AveTime='AveAllR';
handles.timeLineaxCtrl=[];
handles.timeLineaxT=[];

% graphics styles update
colormap(handles.axesCurrent,handles.defColormap);
handles.axesCurrent.CLim=[-10,10];
set(handles.axesCurrent,'XTick',[])
set(handles.axesCurrent,'YTick',[])
set(handles.axesCurrent,'Box','on');

colormap(handles.axesOff,handles.defColormap);
handles.axesOff.CLim=[-10,10];
set(handles.axesOff,'XTick',[]);
set(handles.axesOff,'YTick',[]);
set(handles.axesOff,'Box','on');

colormap(handles.axesSource,handles.defColormap);
handles.axesSource.CLim=[-10,10];
set(handles.axesSource,'XTick',[]);
set(handles.axesSource,'YTick',[]);
set(handles.axesSource,'Box','on');
colorbar(handles.axesSource,'Location','eastoutside');

colormap(handles.axesConverter,copper);
handles.axesConverter.CLim=[0,1];
set(handles.axesConverter,'XTick',[]);
set(handles.axesConverter,'YTick',[]);
set(handles.axesConverter,'Box','on');

colormap(handles.axesConverterLinear,copper);
handles.axesConverterLinear.CLim=[0,1];
set(handles.axesConverterLinear,'XTick',[]);
set(handles.axesConverterLinear,'YTick',[]);
set(handles.axesConverterLinear,'Box','on');

colormap(handles.axesCentroid,copper);
handles.axesCentroid.CLim=[0,1];
set(handles.axesCentroid,'XTick',[]);
set(handles.axesCentroid,'YTick',[]);
set(handles.axesCentroid,'Box','on');

colormap(handles.axesDigInput,copper);
handles.axesDigInput.CLim=[0,1];
set(handles.axesDigInput,'XTick',[]);
set(handles.axesDigInput,'YTick',[]);
set(handles.axesDigInput,'Box','on');

colormap(handles.axesDigOutput,copper);
handles.axesDigOutput.CLim=[0,1];
set(handles.axesDigOutput,'XTick',[]);
set(handles.axesDigOutput,'YTick',[]);
set(handles.axesDigOutput,'Box','on');

colormap(handles.axesTarget,copper);
handles.axesTarget.CLim=[0,1];
set(handles.axesTarget,'XTick',[]);
set(handles.axesTarget,'YTick',[]);
set(handles.axesTarget,'Box','on');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AssocAnalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AssocAnalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in OpenButton.
function OpenButton_Callback(hObject, eventdata, handles)
% hObject    handle to OpenButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.tdms', 'Select an AssocLearningFile');
if isequal(filename,0)
   disp('User selected Cancel')
else
    if ~contains(filename,'AssociativeSwitching')
        disp('Invalid file');
    else
        fpath=fullfile(pathname,filename);
        ref=OpenTFAssoc(fpath);
        if ~isempty(ref)
            %deploy data
            handles.refList={ref};
            handles.PatternList.String=ref.PatternInfo;
            handles.SourceList.String=vertcat(ref.channelsSource,ref.channelsOff);
            handles.ChannelList.String=ref.channels;
            handles.PatternList.Value=1;
            [handles.WaitEdit.String,handles.TrainEdit.String,handles.TestEdit.String]=...
                getWTTText(ref.textdata);
            handles.MetaDataPop.String=fieldnames(ref.textdata);
            handles.SignalEdit.String=ref.textdata.SigFList;
            handles.ControllerEdit.String=ref.textdata.Controllers;
            handles.MarkersEdit.String=ref.textdata.MrkrList;
            handles.ScoringEdit.String=ref.textdata.ScoreList;
            handles.PattPop.String=cellfun(@num2str,num2cell(unique(ref.DtarNumb)),...
                'UniformOutput',false)';
            UpdateSlider(handles);
                
        end
    end
  
end
guidata(hObject,handles);

% --- Executes on selection change in MetaDataPop.
function MetaDataPop_Callback(hObject, eventdata, handles)
% hObject    handle to MetaDataPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 refList=handles.refList{1};
 Sel=handles.MetaDataPop.String{handles.MetaDataPop.Value};
 handles.MetaEdit.String=refList.textdata.(Sel);
 guidata(hObject,handles);
 
% Hints: contents = cellstr(get(hObject,'String')) returns MetaDataPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MetaDataPop

% --- Executes on selection change in PatternList.
function PatternList_Callback(hObject, eventdata, handles)
% hObject    handle to PatternList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Sel=handles.PatternList.Value(1); %only first selected value is taken
UpdateSlider(handles);
handles=plotMultiAxes(handles,Sel,1:9,1:18);

guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns PatternList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PatternList
% --- Executes on slider movement.


function CMapPop_Callback(hObject, eventdata, handles)
% hObject    handle to CMapPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CMapPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CMapPop
handles.defColormap=handles.CMapPop.String{handles.CMapPop.Value};
guidata(hObject,handles);


% --- Executes on button press in RelScalRadio.
function RelScalRadio_Callback(hObject, eventdata, handles)
% hObject    handle to RelScalRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RelScalRadio
tog=get(hObject,'Value');
if isequal(tog,1)
handles.Scaling='scaled';

end
guidata(hObject,handles);

% --- Executes on button press in PlotAllRadio.
function PlotAllRadio_Callback(hObject, eventdata, handles)
% hObject    handle to PlotAllRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotAllRadio
tog=get(hObject,'Value');
if isequal(tog,1)
handles.PlotType='PAll';

end
guidata(hObject,handles);

% --- Executes on button press in PlotSelRadio.
function PlotSelRadio_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSelRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotSelRadio
tog=get(hObject,'Value');
if isequal(tog,1)
handles.PlotType='PSel';

end
guidata(hObject,handles);

% --- Executes on button press in CloneRadio.
function CloneRadio_Callback(hObject, eventdata, handles)
% hObject    handle to CloneRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CloneRadio
tog=get(hObject,'Value');
if isequal(tog,1)
handles.PlotType='PClone';

end
guidata(hObject,handles);

% --- Executes on button press in AbsScalRadio.
function AbsScalRadio_Callback(hObject, eventdata, handles)
% hObject    handle to AbsScalRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AbsScalRadio
tog=get(hObject,'Value');
if isequal(tog,1)
    handles.Scaling='direct';
end
guidata(hObject,handles);

% --- Executes on selection change in ChannelList.
function ChannelList_Callback(hObject, eventdata, handles)
% hObject    handle to ChannelList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChannelList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChannelList
SelChan=hObject.Value;
SelSource=handles.SourceList.Value;
Sel=handles.PatternList.Value;
switch handles.PlotType
    case 'PAll'
        SelChan=1:9;SelSource=1:18;
    case 'PSel'
        
    case 'PClone'
        SelSource=[SelChan SelChan+9];
        handles.SourceList.Value=SelSource;
       
    
end
handles=plotMultiAxes(handles,Sel,SelChan,SelSource);

guidata(hObject,handles);

% --- Executes on button press in UptoR.
function UptoR_Callback(hObject, eventdata, handles)
% hObject    handle to UptoR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UptoR

function TimeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=handles.TimeSlider.Value;
val=round(val);
cT=handles.TimeSlider.UserData(val);
x=[cT cT];
y=handles.axesT.YLim;
y2=handles.axesControl.YLim;

delete([handles.timeLineaxCtrl handles.timeLineaxT]);

handles.timeLineaxCtrl=line(handles.axesControl,x,y2,'Color','r','LineWidth',0.9);
handles.timeLineaxT=line(handles.axesT,x,y,'Color','r','LineWidth',0.9);
guidata(hObject,handles);
ChannelList_Callback(handles.ChannelList,[], handles)

% --- Executes when selected object is changed in FrameToSliderGroup.
function FrameToSliderGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in FrameToSliderGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.AveTime=hObject.Tag;
guidata(hObject,handles);

% --- Executes on button press in SelTestButton.
function SelTestButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelTestButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ref=handles.refList{1};
TrTs=ref.TrTs;
sel=find(cellfun(@(c) isequal(c,'Ts'),TrTs));
handles.PatternList.Value=sel;
guidata(hObject,handles);


% --- Executes on button press in SelTrainButton.
function SelTrainButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelTrainButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ref=handles.refList{1};
TrTs=ref.TrTs;
sel=find(cellfun(@(c) isequal(c,'Tr'),TrTs));
handles.PatternList.Value=sel;
guidata(hObject,handles);

% --- Executes on button press in SelPatternButton.
function SelPatternButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelPatternButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ref=handles.refList{1};
TrTs=ref.TrTs;
TarL=ref.DtarNumb;
CPat=str2double(handles.PattPop.String{handles.PattPop.Value});

SelT=find(TarL==CPat);
SelC=handles.SelPatPop.Value;
switch SelC
    case 1
        sel=find(cellfun(@(c) isequal(c,'Tr'),TrTs));
    case 2
        sel=find(cellfun(@(c) isequal(c,'Ts'),TrTs));
    case 3
        sel=1:length(TrTs);
end

sels=intersect(sel,SelT);
if ~isempty(sels)
handles.PatternList.Value=sels;
end
guidata(hObject,handles);



% --- Executes on selection change in ScoringPop.
function ScoringPop_Callback(hObject, eventdata, handles)
% hObject    handle to ScoringPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ScoringPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ScoringPop


% --- Executes during object creation, after setting all properties.
function ScoringPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScoringPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in ClearAllButton.
function ClearAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ClearSelButton.
function ClearSelButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearSelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function ChannelList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SourceList.
function SourceList_Callback(hObject, eventdata, handles)
% hObject    handle to SourceList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SourceList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SourceList


% --- Executes during object creation, after setting all properties.
function SourceList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SourceList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SelPatPop.
function SelPatPop_Callback(hObject, eventdata, handles)
% hObject    handle to SelPatPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelPatPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelPatPop


% --- Executes during object creation, after setting all properties.
function SelPatPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelPatPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in PattPop.
function PattPop_Callback(hObject, eventdata, handles)
% hObject    handle to PattPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PattPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PattPop


% --- Executes during object creation, after setting all properties.
function PattPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PattPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxKEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MaxKEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxKEdit as text
%        str2double(get(hObject,'String')) returns contents of MaxKEdit as a double


% --- Executes during object creation, after setting all properties.
function MaxKEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxKEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Scorebutton.
function Scorebutton_Callback(hObject, eventdata, handles)
% hObject    handle to Scorebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Sel=handles.PatternList.Value;
SelChan=1:9;
handles=plotMultiScore(handles,Sel,SelChan);
guidata(hObject,handles);



function NumbToConvert_Callback(hObject, eventdata, handles)
% hObject    handle to NumbToConvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumbToConvert as text
%        str2double(get(hObject,'String')) returns contents of NumbToConvert as a double
numb=str2double(hObject.String);
if numb<0
    numb=0;
end

boolarr=NumbToBoolArray(numb);
ColorMapDig=reshape(boolarr,3,3)';
imagesc(handles.axesConverter,ColorMapDig,'CDataMapping','scaled');
caxis(handles.axesConverter,[0 1]);
set(handles.axesConverter,'XTick',[]);
set(handles.axesConverter,'YTick',[]);

imagesc(handles.axesConverterLinear,boolarr,'CDataMapping','scaled');
caxis(handles.axesConverterLinear,[0 1]);
set(handles.axesConverterLinear,'XTick',[]);
set(handles.axesConverterLinear,'YTick',[]);

% --- Executes during object creation, after setting all properties.
function NumbToConvert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumbToConvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
