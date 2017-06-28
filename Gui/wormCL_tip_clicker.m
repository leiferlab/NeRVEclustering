function varargout = wormCL_tip_clicker(varargin)
% WORMCL_TIP_CLICKER MATLAB code for wormCL_tip_clicker.fig
%      WORMCL_TIP_CLICKER, by itself, creates a new WORMCL_TIP_CLICKER or raises the existing
%      singleton*.
%
%      H = WORMCL_TIP_CLICKER returns the handle to a new WORMCL_TIP_CLICKER or the handle to
%      the existing singleton*.
%
%      WORMCL_TIP_CLICKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORMCL_TIP_CLICKER.M with the given input arguments.
%
%      WORMCL_TIP_CLICKER('Property','Value',...) creates a new WORMCL_TIP_CLICKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wormCL_tip_clicker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wormCL_tip_clicker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wormCL_tip_clicker

% Last Modified by GUIDE v2.5 28-Jun-2017 15:54:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @wormCL_tip_clicker_OpeningFcn, ...
    'gui_OutputFcn',  @wormCL_tip_clicker_OutputFcn, ...
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


% --- Executes just before wormCL_tip_clicker is made visible.
function wormCL_tip_clicker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wormCL_tip_clicker (see VARARGIN)

% Choose default command line output for wormCL_tip_clicker
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @show_image);


setappdata(handles.figure1,'holdaxes',false);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'Max',2000);
set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wormCL_tip_clicker wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wormCL_tip_clicker_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function getMousePositionOnImage(src, event,lag)
if nargin==2
    lag=0;
end


handles = guidata(src);
current_axis=handles.axes1;
cursorPoint = get(current_axis, 'CurrentPoint');
curX = cursorPoint(1,1);
curY = cursorPoint(1,2);

xLimits = get(current_axis, 'xlim');
yLimits = get(current_axis, 'ylim');

if (curX > min(xLimits) && curX < max(xLimits) &&...
        curY > min(yLimits) && curY < max(yLimits))
    if handles.get_head.Value
        addPoint('head',curX,curY,handles,lag);
        forward1_Callback(handles.forward1, event, handles)
    elseif handles.get_tail.Value
        addPoint('tail',curX,curY,handles,lag);
        forward1_Callback(handles.forward1, event, handles)
    end
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
    currentData=uipickfiles('FilterSpec',mostRecent,...
        'Prompt','Select the cam1.avi file from the LowMag folder');
catch
    currentData=uipickfiles(...
         'Prompt','Select the cam1.avi file from the LowMag folder');
end
currentData=currentData{1};

setappdata(0,'mostRecent',fileparts(currentData));
set(handles.currentFolder,'String',currentData);

Fid= VideoReader(currentData);
nFrames=Fid.NumberOfFrames;
setappdata(handles.figure1,'Fid',Fid);
setappdata(handles.figure1,'aviFlag',1);


setappdata(handles.figure1,'maxC',0);
setappdata(handles.figure1,'nFrames',nFrames);
set(handles.slider1,'Value',1);
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',nFrames);
set(handles.maxSlider,'String',num2str(nFrames));
set(handles.minSlider,'String','1');
setappdata(handles.figure1,'MaxValue',nFrames);

head_pts=zeros(nFrames,2);
setappdata(handles.figure1,'head_pts',head_pts);
tail_pts=zeros(nFrames,2);
setappdata(handles.figure1,'tail_pts',tail_pts);
show_image(hObject)



function show_image(hObject,eventdata)
%get current image from slider
handles=guidata(get(hObject,'Parent'));
frameNumber=get(handles.slider1,'Value');

%force it to be a multiple of the stepsize
stepSize=str2double(get(handles.stepSize,'String'));
frameNumber=ceil(frameNumber/stepSize)*stepSize;
frameNumber=max(1,round(frameNumber));
frameNumber=min(frameNumber,handles.slider1.Max);
set(handles.slider1,'Value',frameNumber);

set(handles.currentFrame,'String',num2str(frameNumber));

Fid=getappdata(handles.figure1,'Fid');
im_handle=getappdata(handles.figure1,'im_handle');

centerline=getappdata(handles.figure1,'centerline');
CLoffset=getappdata(handles.figure1,'CLoffset');
CLnumber=frameNumber+CLoffset;
if getappdata(handles.figure1,'aviFlag')
    C=read(Fid,frameNumber);
    C=C(:,:,1);
end

maxC=getappdata(handles.figure1,'maxC');
setappdata(handles.figure1,'maxC',max(max(C(:)),maxC));
if isempty(im_handle)
    h=imagesc(C,'Parent',handles.axes1);
    setappdata(handles.figure1,'im_handle',h)
else
    im_handle.CData=C;
end

hold(handles.axes1,'on')
if any(centerline)
    if  get(handles.transpose,'Value')
        X=centerline(:,2,CLnumber);
        Y=centerline(:,1,CLnumber);
    else
        X=centerline(:,1,CLnumber);
        Y=centerline(:,2,CLnumber);
    end
    CL_handle=getappdata(handles.figure1,'CL_handle');
    if ~isempty(CL_handle)
        CL_handle.XData=X;
        CL_handle.YData=Y;
    elseif CLnumber<size(centerline,3)
        CL_handle=plot(handles.axes1,X, Y,'r');
        setappdata(handles.figure1,'CL_handle',CL_handle)
    end    
end

%plot tips
head_handle=getappdata(handles.figure1,'head_handle');

head_pts=getappdata(handles.figure1,'head_pts');
head=head_pts(frameNumber,:);
if ~isempty(head_handle)
    head_handle.XData=head(1);
    head_handle.YData=head(2);
else
    head_handle=plot(handles.axes1,head(1),head(2),'ob');
end

tail_handle=getappdata(handles.figure1,'tail_handle');
tail_pts=getappdata(handles.figure1,'tail_pts');
tail=tail_pts(frameNumber,:);
if ~isempty(tail_handle)
    tail_handle.XData=tail(1);
    tail_handle.YData=tail(2);
else
    tail_handle=plot(handles.axes1,tail(1),tail(2),'og');
end
setappdata(handles.figure1,'head_handle',head_handle);
setappdata(handles.figure1,'tail_handle',tail_handle);
hold(handles.axes1,'off')

drawnow


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
show_image(hObject)


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
show_image(hObject)
%end


% --- Executes on selection change in colorMap.
function colorMap_Callback(hObject, eventdata, handles)
% hObject    handle to colorMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colorMap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorMap
maxC=getappdata(handles.figure1,'maxC');
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
        colormap(handles.axes1,parula(64))
        caxis(handles.axes1,[0,maxC]);
    otherwise
end

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
currentImage=handles.axes1.Children(1).CData;
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
show_image(hObject)

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
elseif  strcmp(eventdata.Key,'f')
    
stepSize=str2double(get(handles.stepSize,'String'));
%user optional delay, if the program is going fast, turn delay on and the
%program will anticipate that you are behind one frame. Works well at max
%speed.

if get(handles.auto_delay,'Value')
    lag=stepSize;
else
    lag=0;
end

    %enter Auto Mode
    display(['Entering Auto Mode, keep your cursor on the head or tail!']);
    display('To exist, press the space bar');
    set(gcf,'Pointer','circle');
    
     while (handles.figure1.CurrentCharacter=='f' ...
             &&  str2double(handles.currentFrame.String)<handles.slider1.Max)
getMousePositionOnImage(handles.axes1, eventdata,lag)
speed=1/get(handles.speed_slider,'Value');
pause(speed)
     end
     
    set(gcf,'Pointer','arrow');

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


% --- Executes on button press in CLselect.
function CLselect_Callback(hObject, eventdata, handles)
% hObject    handle to CLselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');
CLfile=uipickfiles('filterspec',mostRecent);
CLfile=CLfile{1};
centerline=load(CLfile);
CLfieldNames=fieldnames(centerline);
CLfieldIdx=cellfun(@(x) ~isempty(strfind(x,'centerline')),CLfieldNames);
CLoffsetIdx=cellfun(@(x) ~isempty(strfind(x,'off')),CLfieldNames);
if any(CLoffsetIdx)
    CLoffset=centerline.(CLfieldNames{CLoffsetIdx});
else
    CLoffset=0;
end
centerline=centerline.(CLfieldNames{CLfieldIdx});
setappdata(handles.figure1,'CLoffset',CLoffset);
setappdata(handles.figure1,'centerline',centerline);


% --- Executes on button press in transpose.
function transpose_Callback(hObject, eventdata, handles)
% hObject    handle to transpose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transpose


% --- Executes on button press in get_head.
function get_head_Callback(hObject, eventdata, handles)
% hObject    handle to get_head (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    hObject.BackgroundColor=[1 .4 .4];
else
    hObject.BackgroundColor=[0.9400 0.9400 0.9400];
end

if get(handles.get_tail,'Value');
    handles.get_tail.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.get_tail.Value=0;
end


% --- Executes on button press in get_tail.
function get_tail_Callback(hObject, eventdata, handles)
% hObject    handle to get_tail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if get(hObject,'Value')
    hObject.BackgroundColor=[1 .4 .4];
else
    hObject.BackgroundColor=[0.9400 0.9400 0.9400];
end

if get(handles.get_head,'Value');
    handles.get_head.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.get_head.Value=0;
end

function addPoint(location,xselect,yselect,handles,lag)

currentFrame=round(get(handles.slider1,'Value'));
switch location
    case 'head'
        head_pts=getappdata(handles.figure1,'head_pts');
        head_pts(currentFrame-lag,:)=[xselect,yselect];
        setappdata(handles.figure1,'head_pts',head_pts);
    case 'tail'
        tail_pts=getappdata(handles.figure1,'tail_pts');
        tail_pts(currentFrame-lag,:)=[xselect,yselect];
        setappdata(handles.figure1,'tail_pts',tail_pts);
end
set(handles.last_click,'String', num2str(currentFrame));
saveWarning(handles,1);


% --- Executes on button press in save_tips.
function save_tips_Callback(hObject, eventdata, handles)
% hObject    handle to save_tips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
head_pts=getappdata(handles.figure1,'head_pts');
tail_pts=getappdata(handles.figure1,'tail_pts');
currentFolder=get(handles.currentFolder,'String');
currentFolder=fileparts(currentFolder);
fileOutput=[currentFolder filesep 'tip_coodinates'];
save(fileOutput,'head_pts','tail_pts');
saveWarning(handles,0);


function saveWarning(handles,status)
% toggle state of save reminder

if status
    set(handles.save_status,'String','SAVE ME!')
    set(handles.save_status,'BackgroundColor',[.8,0, 0]);
else
    set(handles.save_status,'String','up to date')
    set(handles.save_status,'BackgroundColor',[0,0.8, 0]);
end


function save_status_Callback(hObject, eventdata, handles)
% hObject    handle to save_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_status as text
%        str2double(get(hObject,'String')) returns contents of save_status as a double


% --- Executes during object creation, after setting all properties.
function save_status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_tips.
function load_tips_Callback(hObject, eventdata, handles)
% hObject    handle to load_tips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFolder=get(handles.currentFolder,'String');
currentFolder=fileparts(currentFolder);

tip_file=uipickfiles('Filterspec', currentFolder);
tip_file=tip_file{1};
tip_data=load(tip_file);
head_pts=tip_data.head_pts;
tail_pts=tip_data.tail_pts;

clicks=any([head_pts,tail_pts],2);
last_click=find(clicks,1,'last');
setappdata(handles.figure1,'head_pts',head_pts);
setappdata(handles.figure1,'tail_pts',tail_pts);
set(handles.last_click,'String', num2str(last_click));
saveWarning(handles,0);


% --------------------------------------------------------------------
function adjust_contrast_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to adjust_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%adjust contrast on handles.axes1
handles.axes1.Children(end).CDataMapping='scaled';
h=imcontrast(handles.axes1);
uiwait(h)
%handles.axes1.Children.CDataMapping='direct';



function currentFrame_Callback(hObject, eventdata, handles)
% hObject    handle to currentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentFrame as text
%        str2double(get(hObject,'String')) returns contents of currentFrame as a double

currentFrame=str2double(get(hObject,'String'));
currentFrame=max(currentFrame,handles.slider1.Min);
currentFrame=min(currentFrame,handles.slider1.Max);
handles.slider1.Value=round(currentFrame);


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 


% --- Executes on slider movement.
function speed_slider_Callback(hObject, eventdata, handles)
% hObject    handle to speed_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value= get(hObject,'Value');

if value<1
    value=1;
elseif value>10
    value=10;
end
set(handles.speed_value,'String',num2str(value));
set(handles.speed_slider,'Value',value);

% --- Executes during object creation, after setting all properties.
function speed_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function speed_value_Callback(hObject, eventdata, handles)
% hObject    handle to speed_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speed_value as text
%        str2double(get(hObject,'String')) returns contents of speed_value as a double
value=str2double(get(hObject,'String'));

if value<.1
    value=.1;
elseif value>10
    value=10;
end
set(hObject,'String',num2str(value));
set(handles.speed_slider,'Value',value);

% --- Executes during object creation, after setting all properties.
function speed_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
