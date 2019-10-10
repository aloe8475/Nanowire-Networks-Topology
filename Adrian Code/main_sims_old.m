function varargout = main_sims(varargin)
% Network Simulator v.1
%
% Adrian Diaz Alvarez
% Data Visualization tool for random nanowire networks.
% Mana-NIMS.
%
%
%
%
%
%
%
%
%
% Write here updates to the code:
% -29-10/2018 Add Graph current visualization: 'GCurrent'
% -30-10/2018 Correct video maker (resolution and aspect ratio)
% -01-11/2018 Improve current visualization in graph layout (with custom
% gcurrmap colormap)
% -02-11/2018 Add ElectrodeSeries plot, to plot conductance through
% electrodes in matrix form (visualization for multielectrode correlational
% measurements). ElectrodeNames in simulationsettings are restricted to
% either drain or source type (numbered)
% -06-11 Add Zdenka model (discrete resistance, formation and break
% lengths)
% -
% -
% -

%
% Last Updated 10-26-2018
% Last Modified by GUIDE v2.5 19-Oct-2018 09:37:04

% MAIN_SIMS MATLAB code for main_sims.fig


% Begin initialization code - DO NOT EDIT
%DO NOT EDIT THIS initialization CODE
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_sims_OpeningFcn, ...
    'gui_OutputFcn',  @main_sims_OutputFcn, ...
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


% --- Executes just before main_sims is made visible.
function main_sims_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_sims (see VARARGIN)

%opening functions for simulations, initialization of important structures
%data types and variables
handles.output = hObject;

Networks={};

%
PlotList=table();
PlotList.Names={'TimeSeries';'NetworkLayout';'GraphLayout';'FrequencySeries';'ElectrodeSeries'};
PlotList.Types={...
    {'Current';'CurrentLog';'Voltage';'Conductance';'ConductanceLog';'Resistance';'ResistanceLog';'JunctionWidth';'JunctionWidthLog';'AvJuncWidth'};...
    {'SelEdges';'ShortestPath';'CurrentLayout';'VoltageLayout';'SwitchLayout';'ResistanceLayout'};...
    {'GLayout';'GPath';'GSTree';'GDegree';'GVoltage';'GCurrent';'GWidth';'GBetween';'GMaxflow';'GResistance'};...
    {'PSD'};...
    {'CondSMat';'CurrSMat';'CondDMat';'CurrDMat'}};



handles.PlotList=PlotList;
handles.PlotStatus='first';
handles=FillPlotList(handles);
handles.Networks=Networks;
handles.NetworkSettings=fillNETWORKSETTINGS(2,[],[]);
handles.SimulationSettings=fillSIMULATIONSETTINGS(2,[],[]);
handles.InfoText.String=strsplit('Press Create Network to create new network. Start simulation to start a simulation on the created network'...
    ,'.')';
% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = main_sims_OutputFcn(hObject, eventdata, handles)
% unused at the moment
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function SimList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in NetworkButton.
function NetworkButton_Callback(hObject, ~, handles)
% Creates a network calling the modal windown callNETWORKSETTINGS to fill
% the setting struct with the relevant parameters to construct the network.
% the network is added to the network list.

delete(handles.PlotPanel.Children);
handles.axesPan=axes(handles.PlotPanel);
axes(handles.axesPan);
handles.NetworkSettings=callNETWORKSETTINGS(handles.NetworkSettings);

NewNet=CreateNetwork(handles.InfoText,handles.NetworkSettings);

ParNet=struct();
ParNet.LayOut=NewNet;
ParNet.NetworkSettings=handles.NetworkSettings;
[ParNet.Graph,ParNet.CrossMat.CX,ParNet.CrossMat.CY]...
    =GenerateGraphFromNetwork(NewNet);
ParNet.MetaData=fillNETWORKMETADATA(ParNet);
ParNet.Name=ParNet.MetaData.Name;
ParNet.Simulations={};
%%% add new network to the list
[handles.Networks,handles.NetList.String,handles.NetList.Value]=AddDeleteDataList('add',...
    handles.Networks,ParNet);


guidata(hObject,handles);

% --- Executes on selection change in NetList.
function NetList_Callback(hObject, eventdata, handles)
% Selects a network from network list, retrieve the stored
% simulations for that network and plot the
% nanowire layout and connections in the panel figure.

delete(handles.PlotPanel.Children);
handles.axesPan=axes(handles.PlotPanel);
axes(handles.axesPan);
handles.NodeList.String='NodeList';
handles.NodeList.Value=1;
if isequal(handles.NetList.String,'NetList')
    return;
end
Network=handles.Networks{handles.NetList.Value};
if ~isempty(Network.Simulations)
    handles.SimList.String=cellfun(@(x) x.Name,Network.Simulations,'UniformOutput',false)';
    if handles.SimList.Value>length(handles.SimList.String)
        handles.SimList.Value=length(handles.SimList.String);
    end
else
    handles.SimList.String='SimList';
    handles.SimList.Value=1;
end
PlotNetwork(Network,0);

% --- Executes on button press in SimButton.
function SimButton_Callback(hObject, eventdata, handles)
% Creates a NewSim struct, and adds it to simulation cell array of network

NetList_Callback(handles.NetList, eventdata, handles);
Network=handles.Networks{handles.NetList.Value};
Network.Domains=ConnectedDomains(Network.Graph);
PlotNetwork(Network,'domains');

handles.SimulationSettings=callSIMULATIONSETTINGS(handles.SimulationSettings);
if isempty(handles.SimulationSettings)
    return;
end

NewSim.Settings=handles.SimulationSettings;

[NewSim.Time,NewSim.Electrodes]=getELECTRODES(handles,NewSim.Settings,Network);
if isempty(NewSim.Electrodes)
    guidata(hObject,handles);
    return;
end

NewSim.SelDomain=Network.Domains{NewSim.Electrodes.DomIndex(1)};
NewSim.SelLayout=GetSelectedLayout(Network,NewSim.SelDomain);



if ~isequal(0,handles.ResumeCheck.Value) && ~isequal(handles.SimList.String,'SimList')
    IndexNet=handles.NetList.Value;
    IndexSim=handles.SimList.Value;
    SelSim=handles.Networks{IndexNet}.Simulations{IndexSim};
    NewSim.LastW=SelSim.Data.Wmat{end};%%if you want to resume previous simulation
else
    NewSim.LastW=[];
end

tic;
[NewSim.Data,NewSim.Name]=LaunchSimulation(NewSim);


handles.InfoText.String=strcat('Simulation completed. Elapsed time is',{' '},num2str(toc),{' '},'seconds');
NewSim.SimInfo=GetSimInfo(NewSim.Data);

[handles.NodeList.String,handles.NodeList.UserData]=GetNodeList(NewSim);
handles.NodeList.UserData.Annotate=handles.LabelElectrodeCheck.Value;
handles.NodeList.Value=1;


[Network.Simulations,handles.SimList.String,handles.SimList.Value]=...
    AddDeleteDataList('add',Network.Simulations,NewSim);

handles.Networks{handles.NetList.Value}=Network;

guidata(hObject,handles);


% --- Executes on selection change in SimList.
function SimList_Callback(hObject, eventdata, handles)
% Selects a simulation from the selected network
% updates slider bar and plots the selected subplot
% layout at time 0.

if isequal(handles.SimList.String,'SimList')
    return;
elseif isequal(handles.NetList.String,'NetList')
    handles.SimList.String='SimList';
    handles.SimList.Value=1;
    return;
end

IndexNet=handles.NetList.Value;
IndexSim=handles.SimList.Value;

SelSim=handles.Networks{IndexNet}.Simulations{IndexSim};

t=SelSim.Data.time;
handles.VideoSlider.Max=length(t);handles.VideoSlider.Min=1;
handles.VideoSlider.SliderStep=[1/(length(t)-1) 1/(length(t)-1)];
handles.VideoSlider.Value=1;

handles.PlotSel=SelPlotFromList(handles);
handles.PlotStatus='first';

[handles.NodeList.String,handles.NodeList.UserData]=GetNodeList(SelSim);
handles.NodeList.UserData.Annotate=handles.LabelElectrodeCheck.Value;
%handles.NodeList.Value=1;

handles.InfoText.String=DisplayInfoText(SelSim);

handles=MultiPlot(handles,IndexNet,IndexSim,1);

guidata(hObject,handles);






% --- Executes during object creation, after setting all properties.
function NetList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function VideoSlider_Callback(hObject, eventdata, handles)
% Interpolates the slider position to a time t in the selected simulation
% and plots the selected plots in the figure for the time t.

idx=round(hObject.Value);

if isequal(handles.SimList.String,'SimList')
    return;
end

IndexNet=handles.NetList.Value;
IndexSim=handles.SimList.Value;




NewSel=SelPlotFromList(handles);

if isequal(NewSel,handles.PlotSel)
    handles.PlotStatus='redraw';
else
    handles.PlotStatus='first';
    %you want to change the kind of plot
end

handles.PlotSel=NewSel;
handles.NodeList.UserData.Annotate=handles.LabelElectrodeCheck.Value;

str=handles.InfoText.String;
str{3}=strcat('TimeIndex:',num2str(idx));
handles.InfoText.String=str;


handles=MultiPlot(handles,IndexNet,IndexSim,idx);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function VideoSlider_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in ExtractButton.
function ExtractButton_Callback(hObject, eventdata, handles)
%Extracts current subplot in main_sims panel to a new figure
% use it to copy the figure to save it as a fig or relevant image format

ax=handles.PlotPanel.Children;
fig=figure;

h=copyobj(ax,fig);
set(h,'Units','normalized');



% --- Executes on button press in VideoButton.
function VideoButton_Callback(hObject, eventdata, handles)
% Record a video of the current subplot from time 0 to final time.
% with name, folder and frame rate that you choose.

[Filename,Path]=uiputfile('*.avi','Choose Video File Name');
f=fullfile(Path,Filename);
handles.PlotSel=[];
%create video
V=VideoWriter(f);
V.Quality=100;V.FrameRate=str2double(handles.FrameRateEdit.String);
open(V);
tarfig=figure('Position',[10 10 1024 720]);

max=handles.VideoSlider.Max;

for i=1:max
    handles.VideoSlider.Value=i;
    VideoSlider_Callback(handles.VideoSlider,eventdata,handles);
    frame=getframe(tarfig);
    writeVideo(V,frame);
end
close(V);

delete(tarfig);

%call simlist again so that it creates the proper axes handles again
SimList_Callback(handles.SimList, eventdata, handles);


% --- Executes on button press in DelNetButton.
function DelNetButton_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% Delete a network from network list (with all its simulations)
Ind=handles.NetList.Value;

[handles.Networks,handles.NetList.String,handles.NetList.Value]=AddDeleteDataList('deleteNet',...
    handles.Networks,Ind);
if isequal(handles.NetList.String,'NetList')
    handles.SimList.String='SimList';
end
guidata(hObject,handles);

% --- Executes on button press in DelSimButton.
function DelSimButton_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% Delete the selected simulation for the selected network.
Ind=handles.SimList.Value;
IndNet=handles.NetList.Value;

[handles.Networks{IndNet}.Simulations,handles.SimList.String,handles.SimList.Value]=AddDeleteDataList('deleteSim',...
    handles.Networks{IndNet}.Simulations,Ind);
guidata(hObject,handles);


% --- Executes on button press in any pop button of main_sims
function Pop_Callback(hObject,~,handles)
% Change subplot layout

NewName=hObject.String{hObject.Value};

% all plot list objects are named: 'Plot%i_%iList'
Target=split(hObject.Tag,'p');
TargetList=strcat('Plot',Target{2},'List');

if isequal(NewName,'NoPlot')
    handles.(TargetList).String=NewName;
    handles.(TargetList).Value=1;
    return
end
PlotList=handles.PlotList;


handles.(TargetList).String=PlotList.Types{strcmp(PlotList.Names,NewName)};
handles.(TargetList).Value=1;



guidata(hObject,handles);

% --- Executes on button press in OpenButton.
function OpenButton_Callback(hObject, eventdata, handles)
% Open a previous network saved as a mat file, the network contains
% all the previous simulations (that have not been deleted)

[FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
f=fullfile(PathName,FileName);
s=load(f);
if ~isfield(s,'SelNet')
    return;
end
ParNet=s.SelNet;
[handles.Networks,handles.NetList.String,handles.NetList.Value]=AddDeleteDataList('add',...
    handles.Networks,ParNet);

guidata(hObject,handles);




% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% Save a selected network and its simulations as a mat file
% in the folder that you choose
IndNet=handles.NetList.Value;
SelNet=handles.Networks{IndNet};
Name=SelNet.Name;
Name=strrep(Name,':','_');
Name=strrep(Name,'/','_');
Name=strcat(Name,'.mat');
SelDir=uigetdir();
if ~isequal(SelDir,0)
    name=fullfile(SelDir,Name);
    save(name,'SelNet','-v7.3');
else
    return;
end
%uisave('SelNet',Name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%























































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THESE FUNCTIONS ARE AUTOMATICALLY GENERATED BY MATLAB, I DONT DELTE THEM
% BECAUSE THEY MIGHT BE USED IN FURTHER UPDATES OF THE CODE.




% Hint: get(hObject,'Value') returns toggle state of SimButton
% --- Executes on selection change in Plot1_1List.
function Plot1_1List_Callback(hObject, eventdata, handles)
% hObject    handle to Plot1_1List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot1_1List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot1_1List


% --- Executes during object creation, after setting all properties.
function Plot1_1List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot1_1List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if isequal(hObject.String,'Plot1_1List')
    return;
end



% --- Executes on selection change in Plot2_2List.
function Plot2_2List_Callback(hObject, eventdata, handles)
% hObject    handle to Plot2_2List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot2_2List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot2_2List
if isequal(hObject.String,'Plot2_2List')
    return;
end

% --- Executes during object creation, after setting all properties.
function Plot2_2List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot2_2List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Plot2_1List.
function Plot2_1List_Callback(hObject, eventdata, handles)
% hObject    handle to Plot2_1List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot2_1List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot2_1List
if isequal(hObject.String,'Plot2_1List')
    return;
end


% --- Executes during object creation, after setting all properties.
function Plot2_1List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot2_1List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Plot1_2List.
function Plot1_2List_Callback(hObject, eventdata, handles)
% hObject    handle to Plot1_2List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot1_2List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot1_2List
if isequal(hObject.String,'Plot1_2List')
    return;
end


% --- Executes during object creation, after setting all properties.
function Plot1_2List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot1_2List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NodeList.
function NodeList_Callback(hObject, eventdata, handles)
% hObject    handle to NodeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NodeList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NodeList


% --- Executes during object creation, after setting all properties.
function NodeList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
