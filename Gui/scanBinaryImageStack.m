function varargout = scanBinaryImageStack(varargin)
% SCANBINARYIMAGESTACK MATLAB code for scanBinaryImageStack.fig
%      SCANBINARYIMAGESTACK, by itself, creates a new SCANBINARYIMAGESTACK or raises the existing
%      singleton*.
%
%      H = SCANBINARYIMAGESTACK returns the handle to a new SCANBINARYIMAGESTACK or the handle to
%      the existing singleton*.
%
%      SCANBINARYIMAGESTACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCANBINARYIMAGESTACK.M with the given input arguments.
%
%      SCANBINARYIMAGESTACK('Property','Value',...) creates a new SCANBINARYIMAGESTACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scanBinaryImageStack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scanBinaryImageStack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scanBinaryImageStack

% Last Modified by GUIDE v2.5 08-Dec-2014 12:58:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scanBinaryImageStack_OpeningFcn, ...
                   'gui_OutputFcn',  @scanBinaryImageStack_OutputFcn, ...
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


% --- Executes just before scanBinaryImageStack is made visible.
function scanBinaryImageStack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scanBinaryImageStack (see VARARGIN)

% Choose default command line output for scanBinaryImageStack
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @showImage);


setappdata(handles.figure1,'holdaxes',false);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'Max',2000);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scanBinaryImageStack wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scanBinaryImageStack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in selectFolder.
function selectFolder_Callback(hObject, eventdata, handles)
% hObject    handle to selectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');
Fid=getappdata(handles.figure1,'Fid');
try
    fclose(Fid);
catch
end

try
currentData=uipickfiles('FilterSpec',mostRecent);
catch
    currentData=uipickfiles();
end
currentData=currentData{1};

setappdata(0,'mostRecent',fileparts(currentData));
set(handles.currentFolder,'String',currentData);
if strfind(currentData,'.dat')
Fid=fopen(currentData);
status=fseek(Fid,0,1);
%get image sizes for dat files

[row,col]=getdatdimensions(currentData);
set(handles.imageRows,'String',num2str(row));
set(handles.imageCols,'String',num2str(col));

nFrames=ftell(Fid)/(2*row*col)-1;
setappdata(handles.figure1,'Fid',Fid);
setappdata(handles.figure1,'aviFlag',0);
elseif strfind(currentData,'.avi')
     Fid= VideoReader(currentData);
     nFrames=Fid.NumberOfFrames;
setappdata(handles.figure1,'Fid',Fid);
setappdata(handles.figure1,'aviFlag',1);
end


setappdata(handles.figure1,'maxC',0);
setappdata(handles.figure1,'nFrames',nFrames);
set(handles.slider1,'Value',1);
set(handles.slider1,'Max',nFrames);
set(handles.maxSlider,'String',num2str(nFrames));
set(handles.minSlider,'String','1');
setappdata(handles.figure1,'MaxValue',nFrames);
showImage(hObject)


function showImage(hObject,eventdata)

handles=guidata(get(hObject,'Parent'));
frameNumber=get(handles.slider1,'Value');
frameNumber=max(1,round(frameNumber));
Fid=getappdata(handles.figure1,'Fid');

if getappdata(handles.figure1,'aviFlag')
    C=read(Fid,frameNumber);
    C=C(:,:,1);
else
    
row=str2double(get(handles.imageRows,'String'));
col=str2double(get(handles.imageCols,'String'));
nPix=row*col;


status=fseek(Fid,2*frameNumber*nPix,-1);
if ~status
  frewind(Fid) 
  status=fseek(Fid,2*frameNumber*nPix,-1);
end
pixelValues=fread(Fid,nPix,'uint16',0,'l');
C=reshape(pixelValues,row,col);
end

maxC=getappdata(handles.figure1,'maxC');
setappdata(handles.figure1,'maxC',max(max(C(:)),maxC));

h=imagesc(C,'Parent',handles.axes1);
set(handles.currentFrame,'String',num2str(frameNumber));
switch get(handles.colorMap,'Value');
    case 1
        colormap(handles.axes1,jet(64))
      %  caxis(handles.axes1,[0,maxC]); 

    case 2
        colormap(handles.axes1,hot(64))
        caxis(handles.axes1,[0,maxC]); 

    case 3
        colormap(handles.axes1,gray(64))
        caxis(handles.axes1,[0,maxC]); 

    case 4
        C(1:1:row/2,:)=normalizeRange((C(1:1:row/2,:)))+1;
        C(round(row/2):end,:)=normalizeRange((C(round(row/2):end,:)));
        h=imagesc(C,'Parent',handles.axes1);
      colormap(handles.axes1,[hot(32);circshift(hot(32),[0,1])]);
      caxis(handles.axes1,[0,2]); 
        
    otherwise
end
drawnow
setappdata(handles.figure1,'currentImage',C);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function minSlider_Callback(hObject, eventdata, handles)
% hObject    handle to minSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSlider as text
%        str2double(get(hObject,'String')) returns contents of minSlider as a double
val=str2double(get(hObject,'String'));
set(handles.slider1,'Min',val);
set(handles.slider1,'Value',val);


% --- Executes during object creation, after setting all properties.
function minSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to maxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxSlider as text
%        str2double(get(hObject,'String')) returns contents of maxSlider as a double
val=str2double(get(hObject,'String'));
set(handles.slider1,'Max',val);

% --- Executes during object creation, after setting all properties.
function maxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in back1.
function back1_Callback(hObject, eventdata, handles)
% hObject    handle to back1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFrame=get(handles.slider1,'Value');
stepSize=str2double(get(handles.stepSize,'String'));
set(handles.slider1,'Value',max(currentFrame-stepSize,1));
showImage(hObject)


% --- Executes on button press in forward1.
function forward1_Callback(hObject, eventdata, handles)
% hObject    handle to forward1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get(handles.forward1)
%while strcmp(get(handles.forward1,'Selected'),'on')
currentFrame=get(handles.slider1,'Value');
stepSize=str2double(get(handles.stepSize,'String'));
set(handles.slider1,'Value',min(currentFrame+stepSize,get(handles.slider1,'max')));
showImage(hObject)
%end



function imageRows_Callback(hObject, eventdata, handles)
% hObject    handle to imageRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageRows as text
%        str2double(get(hObject,'String')) returns contents of imageRows as a double


% --- Executes during object creation, after setting all properties.
function imageRows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imageCols_Callback(hObject, eventdata, handles)
% hObject    handle to imageCols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageCols as text
%        str2double(get(hObject,'String')) returns contents of imageCols as a double


% --- Executes during object creation, after setting all properties.
function imageCols_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageCols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colorMap.
function colorMap_Callback(hObject, eventdata, handles)
% hObject    handle to colorMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colorMap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorMap


% --- Executes during object creation, after setting all properties.
function colorMap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function snapShotName_Callback(hObject, eventdata, handles)
% hObject    handle to snapShotName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snapShotName as text
%        str2double(get(hObject,'String')) returns contents of snapShotName as a double
setappdata(handles.figure1,'ender',1);

% --- Executes during object creation, after setting all properties.
function snapShotName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snapShotName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in snapshot.
function snapshot_Callback(hObject, eventdata, handles)
% hObject    handle to snapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentImage=getappdata(handles.figure1,'currentImage');
currentFolder=getappdata(0,'mostRecent');
imageName=get(handles.snapShotName,'String');
imageName=fullfile(currentFolder,imageName);
ender=getappdata(handles.figure1,'ender');

imageName=[imageName,num2str(ender,'%3.5d') '.tif'];


if isempty(ender);
    ender=1;
else
    ender=ender+1;
end

setappdata(handles.figure1,'ender',ender);

tiffwrite(imageName,single(currentImage),'tif',0);


currentFrame=get(handles.slider1,'Value');
set(handles.slider1,'Value',min(currentFrame+1,get(handles.slider1,'max')));
showImage(hObject)

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


        if strcmp(eventdata.Key,'rightarrow')||strcmp(eventdata.Key,'d')%strcmp(evnt.Key,'space') || 
            forward1_Callback(handles.slider1,eventdata,handles);
            
        %Backward
        elseif strcmp(eventdata.Key,'leftarrow')||strcmp(eventdata.Key,'a')
            back1_Callback(handles.slider1,eventdata,handles);
        %Up
        
        elseif  strcmp(eventdata.Key,'space')
            snapshot_Callback(handles.slider1,eventdata,handles);
             forward1_Callback(handles.slider1,eventdata,handles);
   
        end



function stepSize_Callback(hObject, eventdata, handles)
% hObject    handle to stepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepSize as text
%        str2double(get(hObject,'String')) returns contents of stepSize as a double


% --- Executes during object creation, after setting all properties.
function stepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function V=normalizeRange(X)
%The function takes a scalar or nd matrix X and normalize the matrix to a
%minimum of 0 and a maximum of 1. If all the values in the matrix are the
%same, they are all replaced with ones. 
if ~isempty(X)
L=numel(X);
XX=reshape(X,1,L);
if max(XX)==min(XX)
    V=ones(size(X));
else
V=(X-min(XX))/(max(XX)-min(XX));
end
else
    V=X;
end