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

% Last Modified by GUIDE v2.5 03-Nov-2016 11:28:59

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
row=str2double(get(handles.imageRows,'String'));
col=str2double(get(handles.imageCols,'String'));

setappdata(0,'mostRecent',fileparts(currentData));
set(handles.currentFolder,'String',currentData);
if strfind(currentData,'.dat')
Fid=fopen(currentData);
status=fseek(Fid,0,1);
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

head_pts=zeros(nFrames,2);
setappdata(handles.figure1,'head_pts',head_pts);
tail_pts=zeros(nFrames,2);
setappdata(handles.figure1,'tail_pts',tail_pts);
show_image(hObject)



function show_image(hObject,eventdata)

handles=guidata(get(hObject,'Parent'));
frameNumber=get(handles.slider1,'Value');
frameNumber=max(1,round(frameNumber));
Fid=getappdata(handles.figure1,'Fid');
centerline=getappdata(handles.figure1,'centerline');
CLoffset=getappdata(handles.figure1,'CLoffset');
CLnumber=frameNumber+CLoffset;
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
if 0
C=pedistalSubtract(C);
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
        C(1:1:row/2,:)=normalizeRange(pedistalSubtract(C(1:1:row/2,:)))+1;
        C(round(row/2):end,:)=normalizeRange(pedistalSubtract(C(round(row/2):end,:)));
        h=imagesc(C,'Parent',handles.axes1);
      colormap(handles.axes1,[hot(32);circshift(hot(32),[0,1])]);
      caxis(handles.axes1,[0,2]); 
    case 5
        colormap(handles.axes1,parula(64))
        caxis(handles.axes1,[0,maxC]); 
        
    otherwise
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
if CLnumber<size(centerline,3)
plot(handles.axes1,X, Y,'r');
end
plot(handles.axes1,X([1,end],:), Y([1 end],1),'or');

end

%plot tips

head_pts=getappdata(handles.figure1,'head_pts');
head=head_pts(frameNumber,:);
if any(head)
    plot(handles.axes1,head(1),head(2),'ob');
end

tail_pts=getappdata(handles.figure1,'tail_pts');
tail=tail_pts(frameNumber,:);
if any(tail)
    plot(handles.axes1,tail(1),tail(2),'og');
end
hold(handles.axes1,'off')

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
        
        elseif strcmp(eventdata.Key,'w')
            get_head_Callback(handles.slider1,eventdata,handles);            
        elseif strcmp(eventdata.Key,'s')
            get_tail_Callback(handles.slider1,eventdata,handles);            
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
currentFrame=get(handles.slider1,'Value');
head_pts=getappdata(handles.figure1,'head_pts');
[xselect,yselect]=ginput(1);
head_pts(currentFrame,:)=[xselect,yselect];
setappdata(handles.figure1,'head_pts',head_pts);
set(handles.last_click,'String', num2str(currentFrame));
saveWarning(handles,1);



% --- Executes on button press in get_tail.
function get_tail_Callback(hObject, eventdata, handles)
% hObject    handle to get_tail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFrame=get(handles.slider1,'Value');
tail_pts=getappdata(handles.figure1,'tail_pts');
[xselect,yselect]=ginput(1);
tail_pts(currentFrame,:)=[xselect,yselect];
setappdata(handles.figure1,'tail_pts',tail_pts);
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
imcontrast(handles.axes1);
