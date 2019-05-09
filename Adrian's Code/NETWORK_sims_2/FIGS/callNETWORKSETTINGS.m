function varargout = callNETWORKSETTINGS(varargin)
% CALLNETWORKSETTINGS MATLAB code for callNETWORKSETTINGS.fig
%      CALLNETWORKSETTINGS, by itself, creates a new CALLNETWORKSETTINGS or raises the existing
%      singleton*.
%
%      H = CALLNETWORKSETTINGS returns the handle to a new CALLNETWORKSETTINGS or the handle to
%      the existing singleton*.
%
%      CALLNETWORKSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALLNETWORKSETTINGS.M with the given input arguments.
%
%      CALLNETWORKSETTINGS('Property','Value',...) creates a new CALLNETWORKSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before callNETWORKSETTINGS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to callNETWORKSETTINGS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help callNETWORKSETTINGS

% Last Modified by GUIDE v2.5 16-Oct-2018 11:43:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @callNETWORKSETTINGS_OpeningFcn, ...
                   'gui_OutputFcn',  @callNETWORKSETTINGS_OutputFcn, ...
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


% --- Executes just before callNETWORKSETTINGS is made visible.
function callNETWORKSETTINGS_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to callNETWORKSETTINGS (see VARARGIN)

% Choose default command line output for callNETWORKSETTINGS
handles.output = hObject;

[handles]=fillNETWORKSETTINGS(0,varargin{1},handles);

% Update handles structure
guidata(hObject, handles);


set(handles.figure1,'WindowStyle','Modal');

% UIWAIT makes callNETWORKSETTINGS wait for user response (see UIRESUME)

 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = callNETWORKSETTINGS_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


varargout{1} = handles.NetworkSettings;
delete(handles.figure1);

% --- Executes on button press in CreateButton.
function CreateButton_Callback(hObject, ~, handles)
% hObject    handle to CreateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[handles.NetworkSettings]=fillNETWORKSETTINGS(1,[],handles);

guidata(hObject,handles);

uiresume(handles.figure1);
