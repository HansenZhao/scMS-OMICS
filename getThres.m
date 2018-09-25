function varargout = getThres(varargin)
% GETTHRES MATLAB code for getThres.fig
%      GETTHRES, by itself, creates a new GETTHRES or raises the existing
%      singleton*.
%
%      H = GETTHRES returns the handle to a new GETTHRES or the handle to
%      the existing singleton*.
%
%      GETTHRES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETTHRES.M with the given input arguments.
%
%      GETTHRES('Property','Value',...) creates a new GETTHRES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getThres_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getThres_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getThres

% Last Modified by GUIDE v2.5 25-Sep-2018 16:32:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getThres_OpeningFcn, ...
                   'gui_OutputFcn',  @getThres_OutputFcn, ...
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


% --- Executes just before getThres is made visible.
function getThres_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getThres (see VARARGIN)

% Choose default command line output for getThres

handles.vec = varargin{1};
handles.mz = varargin{2};
handles.mat = varargin{3};
handles.slider.Min = min(handles.vec);
handles.slider.Max = max(handles.vec);
handles.slider.Value = handles.slider.Min;
handles.curValue = handles.slider.Min;
refreshGUI(handles);
% Update handles structure
guidata(hObject, handles);
handles.output = handles.curValue;
% UIWAIT makes getThres wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = getThres_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.curValue;
close(handles.figure1);


% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curValue = get(hObject,'Value');
refreshGUI(handles);
guidata(hObject,handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edt_thres_Callback(hObject, eventdata, handles)
% hObject    handle to edt_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2double(get(hObject,'String'));
if tmp>=handles.slider.Min && tmp<=handles.slider.Max
    handles.curValue = tmp;
    handles.slider.Value = tmp;
else
    hObject.String = num2str(handles.curValue);
end
guidata(hObject,handles);
refreshGUI(handles);
% Hints: get(hObject,'String') returns contents of edt_thres as text
%        str2double(get(hObject,'String')) returns contents of edt_thres as a double


% --- Executes during object creation, after setting all properties.
function edt_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_OK.
function btn_OK_Callback(hObject, eventdata, handles)
% hObject    handle to btn_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

% --- Executes on button press in Btn_Cancel.
function Btn_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Btn_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curValue = -1;
guidata(hObject,handles);
uiresume(handles.figure1);

function refreshGUI(handles)
    handles.edt_thres.String = num2str(handles.slider.Value);
    I = handles.vec > handles.curValue;
    y = mean(handles.mat(:,I),1);
    x = handles.mz(I);
    drawMZWithTag(handles.axes1,x,y,20);
    plot(handles.axes2,handles.mz,handles.vec,'LineWidth',1);
    handles.txt_begin.String = num2str(handles.slider.Min);
    handles.txt_end.String = num2str(handles.slider.Max);
    handles.txt_cur.String = num2str(handles.slider.Value);
    handles.axes1.XLim = [min(handles.mz),max(handles.mz)];
    handles.axes2.XLim = [min(handles.mz),max(handles.mz)];
    handles.txt_num.String = num2str(sum(I));
