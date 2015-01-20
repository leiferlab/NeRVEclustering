function varargout = VisualizeWormDataFiducialAligned(varargin)
% VISUALIZEWORMDATAFIDUCIALALIGNED MATLAB code for VisualizeWormDataFiducialAligned.fig
%      VISUALIZEWORMDATAFIDUCIALALIGNED, by itself, creates a new VISUALIZEWORMDATAFIDUCIALALIGNED or raises the existing
%      singleton*.
%
%      H = VISUALIZEWORMDATAFIDUCIALALIGNED returns the handle to a new VISUALIZEWORMDATAFIDUCIALALIGNED or the handle to
%      the existing singleton*.
%
%      VISUALIZEWORMDATAFIDUCIALALIGNED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZEWORMDATAFIDUCIALALIGNED.M with the given input arguments.
%
%      VISUALIZEWORMDATAFIDUCIALALIGNED('Property','Value',...) creates a new VISUALIZEWORMDATAFIDUCIALALIGNED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisualizeWormDataFiducialAligned_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisualizeWormDataFiducialAligned_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisualizeWormDataFiducialAligned

% Last Modified by GUIDE v2.5 23-Dec-2014 16:06:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VisualizeWormDataFiducialAligned_OpeningFcn, ...
    'gui_OutputFcn',  @VisualizeWormDataFiducialAligned_OutputFcn, ...
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


% --- Executes just before VisualizeWormDataFiducialAligned is made visible.
function VisualizeWormDataFiducialAligned_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualizeWormDataFiducialAligned (see VARARGIN)

% Choose default command line output for VisualizeWormDataFiducialAligned

%set up slider
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @plotter);
hlistenerz=addlistener(handles.zSlider,'ContinuousValueChange',...
    @plotter);
setappdata(handles.zSlider,'hlistener',hlistenerz);
set(handles.zSlider,'SliderStep',[1,1]);

setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'SliderStep',[1,1]);

hlistener2=addlistener(handles.slider2,'ContinuousValueChange',...
    @plotSlide);
setappdata(handles.slider2,'hlistener',hlistener2);
set(handles.slider2,'SliderStep',[1,1]);



% set up timer for video play
playt.TimerFcn = {@TmrFcn,handles};
playt.BusyMode = 'Queue';
playt.ExecutionMode = 'FixedRate';
playt.Period = 1/2; % set this to 2, then will make the play skip frames accordingly
setappdata(handles.figure1,'playt',playt);
setappdata(handles.figure1,'playTimer',timer(playt));
setappdata(handles.figure1,'FPS',1)
set(handles.framesPerSecond,'String','1')



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VisualizeWormDataFiducialAligned wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualizeWormDataFiducialAligned_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SelectFolder.
function SelectFolder_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%select folers with data,

%get last loaded folder if there is one
% imFolder=getappdata(0,'imFolder');
% 
% %load mat file folder
% if isempty(imFolder)
%     imFolder = uipickfiles('prompt','Select image Folders');
%     setappdata(0,'imFolder',imFolder{1});
% else
%     try
%     imFolder = uipickfiles('filterspec',imFolder,'prompt','Select image Folders');
%         setappdata(0,'imFolder',imFolder{1});
%     catch
%     imFolder = uipickfiles('prompt','Select image Folders');
%         setappdata(0,'imFolder',imFolder{1});
%     end
% end

%select image folder
display('Select image folder, select one folder for split image or select red then green folder');
rawImFolder=uipickfiles;%('FilterSpec', fileparts(imFolder));


if isdir(rawImFolder{1})
for iFolder=1:length(rawImFolder)
imFiles{iFolder}=dir([rawImFolder{iFolder} filesep '*.tif']);
end
end

firstFile=imFiles{1};
stackName=firstFile(1).name;
imageInfo=imfinfo([rawImFolder{1} filesep stackName]);

setappdata(0,'imFiles',imFiles);
setappdata(0,'rawImFolder',rawImFolder);

%setting slider parameters
set(handles.slider1,'Min',1)
if isempty(imFiles)
    set(handles.slider1,'Max',2);
else
    set(handles.slider1,'Max',length(imFiles{1}));
end

set(handles.slider1,'Value',1)


set(handles.maxTime,'String',length(imFiles{1}));
set(handles.minTime,'String',1);

setappdata(handles.figure1,'currentFrame',1);
set(handles.slider1,'value',1);
%set(handles.slider2,'min',1);
%set(handles.slider2,'max',max(trackOutput(:,end)));
%set(handles.slider2,'value',1);
set(handles.zSlider,'value',1);
set(handles.zSlider,'min',1);
set(handles.zSlider,'max',length(imageInfo));

plotter(handles.slider1,eventdata);








% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%plotter(handles.slider1,eventdata);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.figure1,'currentFrame',get(handles.slider1,'value'));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in channelSelect.
function channelSelect_Callback(hObject, eventdata, handles)
% hObject    handle to channelSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,eventdata);
% Hints: contents = cellstr(get(hObject,'String')) returns channelSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channelSelect


% --- Executes during object creation, after setting all properties.
function channelSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channelSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function plotter(hObject,eventdata)

%plots current frame in image window and plot window

handles=guidata(get(hObject,'Parent'));
timeStep=str2double(get(handles.timeStep,'string'));
smoothWindow=str2double(get(handles.smoothingWindow,'String'));
startTime=str2double(get(handles.startTime,'String'));
normalizeFlag=get(handles.normalizeButton,'value');
imFolder=getappdata(0,'imFolder');
iImage=round(get(handles.slider1,'Value'));
minFrame=str2double(get(handles.minTime,'string'));
maxFrame=str2double(get(handles.maxTime,'string'));
minFrame=round(minFrame/timeStep);
maxFrame=round(maxFrame/timeStep);
zSlice=round(get(handles.zSlider,'value'));

%load track data, and lists of mat and tif file
imFiles=getappdata(0,'imFiles');
rawImFolder=getappdata(0,'rawImFolder');
                stackName=imFiles{1};
                stackName=stackName(iImage).name;

        switch get(handles.channelSelect,'value')
            case 1
                baseImg= pedistalSubtract(double(imread([rawImFolder{1} filesep stackName]...
                    ,'tif','index',zSlice))); %red
            case 2
                baseImg= pedistalSubtract(double(imread([rawImFolder{2} filesep stackName]...
                    ,'tif','index',zSlice))); %green
            case 3
                baseImg=labelMask;
            case 4
                baseImg=0;
                for iSlice=1:length(zPos)/2
                    imSlice=double(imread([rawImFolder{1} filesep stackName]...
                        ,'tif','index',iSlice));
                    imSlice=pedistalSubtract(imSlice);
                    baseImg=nanmax( baseImg,imSlice);
                end
                baseImg=baseImg;
            otherwise
                baseImg=[];
           %     centroids(:,3)=centroids(:,3)-mean(centroids(:,3))+50;
        end
        

setappdata(0,'baseImg',baseImg)
set(handles.FrameIdx,'String', [ stackName '   ' num2str(zSlice)]);
imagesc(baseImg,'parent',handles.axes1')

        newContrast=getappdata(handles.figure1,'newContrast');
        if isempty(newContrast)
            newContrast=[min(baseImg(:)),max(baseImg(:))];
        end
%         baseImg(baseImg<newContrast(1)) = newContrast(1);
%         baseImg(baseImg>newContrast(2)) = newContrast(2);
%         baseImg = (baseImg-newContrast(1))./diff(newContrast);
        hold(handles.axes1,'off')

caxis(handles.axes1, [newContrast]);
%     figure



function smoothingWindow_Callback(hObject, eventdata, handles)
% hObject    handle to smoothingWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,eventdata);

% Hints: get(hObject,'String') returns contents of smoothingWindow as text
%        str2double(get(hObject,'String')) returns contents of smoothingWindow as a double


% --- Executes during object creation, after setting all properties.
function smoothingWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothingWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function startTime_Callback(hObject, eventdata, handles)
% hObject    handle to startTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,eventdata);
% Hints: get(hObject,'String') returns contents of startTime as text
%        str2double(get(hObject,'String')) returns contents of startTime as a double


% --- Executes during object creation, after setting all properties.
function startTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cmapping.
function cmapping_Callback(hObject, eventdata, handles)
% hObject    handle to cmapping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
Cstyle=contents{get(hObject,'Value')};
colormap(handles.axes1,Cstyle);

% Hints: contents = cellstr(get(hObject,'String')) returns cmapping contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmapping


% --- Executes during object creation, after setting all properties.
function cmapping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmapping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in goBack.
function goBack_Callback(hObject, eventdata, handles)
% hObject    handle to goBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FPS=getappdata(handles.figure1,'FPS');
set(handles.slider1,'value',get(handles.slider1,'value')-FPS);
setappdata(handles.figure1,'currentFrame',get(handles.slider1,'value'));
plotter(handles.slider1,eventdata)


% --- Executes on button press in goForward.
function goForward_Callback(hObject, eventdata, handles)
% hObject    handle to goForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FPS=getappdata(handles.figure1,'FPS');
set(handles.slider1,'value',get(handles.slider1,'value')+FPS);
setappdata(handles.figure1,'currentFrame',get(handles.slider1,'value'));
plotter(handles.slider1,eventdata)




% --- Executes on button press in playVideo.
function playVideo_Callback(hObject, eventdata, handles)
% hObject    handle to playVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
% disp(button_state);
% disp(get(hObject,'Max'));
% disp(get(hObject,'Min'));
if button_state == get(hObject,'Max')
    % Toggle button is pressed, take appropriate action
    set(hObject,'String','Stop');
    set(hObject,'ForegroundColor',[1 0 0]);
    
    %     set(handles.cursortoggle,'State','off'); % having cursor on creates errs
    start(getappdata(handles.figure1,'playTimer'))
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed, take appropriate action
    set(hObject,'String','Play');
    set(hObject,'ForegroundColor',[0 1 0]);
    
    stop(getappdata(handles.figure1,'playTimer'))
    %     set(handles.cursortoggle,'State','on'); % cursor on is default!
end

% Hint: get(hObject,'Value') returns toggle state of playVideo



function framesPerSecond_Callback(hObject, eventdata, handles)
% hObject    handle to framesPerSecond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FPS=str2double(get(handles.framesPerSecond,'String'));
setappdata(handles.figure1,'FPS',FPS)


% Hints: get(hObject,'String') returns contents of framesPerSecond as text
%        str2double(get(hObject,'String')) returns contents of framesPerSecond as a double


% --- Executes during object creation, after setting all properties.
function framesPerSecond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framesPerSecond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TmrFcn(src,event,handles)
% pull appdata from the handles structure

CurrentFrame = getappdata(handles.figure1,'currentFrame');
% set(handles.slider1,'Value',CurrentFrame)
%CurrentFrame=1;
totalFrames = get(handles.slider1,'Max');
frameStep=str2double(get(handles.framesPerSecond,'string'))/2;
setappdata(handles.figure1,'currentFrame',CurrentFrame+frameStep);
loop = false;

% at some point, include the ability to loop
% loop = get(handles.loop,'Value');
% if the current frame is less than the total, increment frame by one
if CurrentFrame < totalFrames
    set(handles.slider1,'Value',floor(CurrentFrame+frameStep))
    % otherwise, if looping, reset to 1
    % elseif loop == get(handles.loop,'Max')
elseif loop
    set(handles.slider1,'Value',1)
    % otherwise, stop playback
else
    set(handles.playVideo_Callback,'Value',get(handles.togglePlay,'Min'));
    playVideo_Callback(handles.togglePlay, event, handles);
end
plotter(handles.slider1,event)
if get(handles.makeMovie,'value')
    frame=getframe(handles.axes1);
    writerObj=getappdata(handles.figure1,'writerObj');
    writeVideo(writerObj,frame)
    
    
    frame2=getframe(handles.axes2);
    writerObj2=getappdata(handles.figure1,'writerObj2');
    writeVideo(writerObj2,frame2)
    
end



% --- Executes on selection change in plotChannel.
function plotChannel_Callback(hObject, eventdata, handles)
% hObject    handle to plotChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,'eventdata');
% Hints: contents = cellstr(get(hObject,'String')) returns plotChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotChannel


% --- Executes during object creation, after setting all properties.
function plotChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectNeuron1.
function selectNeuron1_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',1);
cursorNeuronSelect(hObject,eventdata)
displayIdx=get(handles.DisplayIdx,'data');
displayIdx{1,1}=round(get(handles.slider2,'value'));
set(handles.slider2,'value',displayIdx{1,1});


% --- Executes on button press in selectNeuron2.
function selectNeuron2_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',2);
cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron3.
function selectNeuron3_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',3);
cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron4.
function selectNeuron4_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',4);
cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron5.
function selectNeuron5_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',5);
cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron6.
function selectNeuron6_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',6);
cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron7.
function selectNeuron7_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',7);
cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron8.
function selectNeuron8_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',8);
cursorNeuronSelect(hObject,eventdata)

function cursorNeuronSelect(hObject,eventdata)
handles=guidata(get(hObject,'Parent'));
plotIdx=getappdata(handles.figure1,'cursorTarget');
currentCentroids=getappdata(handles.figure1,'currentCentroids');

[xselect,yselect]=ginput(1);
xRange=xlim(handles.axes1);
yRange=ylim(handles.axes1);
if xselect>xRange(1) && xselect< xRange(2) && yselect>yRange(1) && yselect<yRange(2);
    minD=pdist2([xselect,yselect],currentCentroids(:,1:2),'euclidean','smallest',1);
    pointIdx=find(minD==min(minD),1,'first');
    pointIdx=currentCentroids(pointIdx,3);
else
    pointIdx=nan;
end

displayIdx=get(handles.DisplayIdx,'data');
displayIdx{plotIdx,1}=pointIdx;
set(handles.DisplayIdx,'data',displayIdx);
plotter(handles.slider1,'eventdata');


% --- Executes on button press in clearPlots.
function clearPlots_Callback(hObject, eventdata, handles)
% hObject    handle to clearPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.DisplayIdx,'data', {[];[];[];[];[];[];[];[]});
cla(handles.axes2);



% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in sortPlots.
function sortPlots_Callback(hObject, eventdata, handles)
% hObject    handle to sortPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
displayIdx=get(handles.DisplayIdx,'data');
plotIdx=[displayIdx{:,1}];
plotIdx=sort(unique(plotIdx));
clearPlots_Callback(hObject,eventdata,handles);
for i =1:length(plotIdx);
    displayIdx{i,1}=plotIdx(i);
end
for i=length(plotIdx)+1:length(displayIdx)
    displayIdx{i,1}=[];
end


set(handles.DisplayIdx,'data',displayIdx);
plotter(handles.slider1,'eventdata');


function plotSlide(hObject,eventdata)
handles=guidata(get(hObject,'Parent'));
displayIdx=get(handles.DisplayIdx,'data');
displayIdx{1,1}=round(get(handles.slider2,'value'));
set(handles.DisplayIdx,'data',displayIdx);
plotter(handles.slider1,'eventdata');


% --- Executes on button press in goForward2.
function goForward2_Callback(hObject, eventdata, handles)
% hObject    handle to goForward2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)set(handles.slider1,'value',get(handles.slider1,'value')+1);
set(handles.slider2,'value',min(get(handles.slider2,'value')+1,get(handles.slider2,'max')));
plotSlide(hObject,eventdata)


% --- Executes on button press in runTrack.
function runTrack_Callback(hObject, eventdata, handles)
% hObject    handle to runTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imFolder=getappdata(0,'imFolder');
trackFiles=dir([imFolder filesep 'track*.mat']);

if numel(trackFiles)~=1
    trackFiles=uipickfiles('filterSpec',imFolder);
    trackFiles=trackFiles{1};
else
    
    trackFiles=[imFolder filesep trackFiles(1).name];
end

load(trackFiles);

setappdata(handles.figure1,'trackOutput',trackOutput);





% --- Executes on selection change in plotChannel2.
function plotChannel2_Callback(hObject, eventdata, handles)
% hObject    handle to plotChannel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,'eventdata');

% Hints: contents = cellstr(get(hObject,'String')) returns plotChannel2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotChannel2


% --- Executes during object creation, after setting all properties.
function plotChannel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotChannel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function movieTitle_Callback(hObject, eventdata, handles)
% hObject    handle to movieTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movieTitle as text
%        str2double(get(hObject,'String')) returns contents of movieTitle as a double


% --- Executes during object creation, after setting all properties.
function movieTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movieTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makeMovie.
function makeMovie_Callback(hObject, eventdata, handles)
% hObject    handle to makeMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state=get(hObject,'value');
ax1movie=[get(handles.movieTitle,'String') '1'];
ax2movie=[get(handles.movieTitle,'String') '2'];
startTime=str2double(get(handles.startTime,'String'));
imFolder=getappdata(0,'imFolder');

if button_state
    writerObj=VideoWriter([imFolder filesep ax1movie]);
    writerObj2=VideoWriter([imFolder filesep ax2movie]);
    
    setappdata(handles.figure1,'writerObj',writerObj);
    setappdata(handles.figure1,'writerObj2',writerObj2);
    
    open(writerObj);
    open(writerObj2);
    
    set(handles.slider1,'value',startTime)
    set(handles.playVideo,'value',1);
    %playVideo_Callback(handles.playVideo, eventdata, handles)
    set(hObject,'String','Rec');
else
    writerObj=getappdata(handles.figure1,'writerObj');
    writerObj2=getappdata(handles.figure1,'writerObj2');
    
    close(writerObj);
    close(writerObj2);
    
    set(handles.playVideo,'value',0);
    playVideo_Callback(handles.playVideo, eventdata, handles)
    set(hObject,'String','Make Movie')
end


% --- Executes on button press in normalizeButton.
function normalizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to normalizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalizeButton
plotter(handles.slider1,eventdata);



function timeStep_Callback(hObject, eventdata, handles)
% hObject    handle to timeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeStep as text
%        str2double(get(hObject,'String')) returns contents of timeStep as a double
timeStep=str2double(get(handles.timeStep,'string'));
maxFrame=str2double(get(handles.maxTime,'String'))/timeStep;
set(handles.maxTime,'String',num2str(maxFrame));
minFrame=str2double(get(handles.minTime,'String'))/timeStep;
set(handles.minTime,'String',num2str(minFrame));


% --- Executes during object creation, after setting all properties.
function timeStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function adjustContrast_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to adjustContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% display the base image, calculate the display image;
baseImg = getappdata(handles.figure1,'baseImg');
imageHandle = findobj('Parent',handles.axes1,'type','image');
storedImage = get(imageHandle);
set(imageHandle,'cdata',baseImg,'cdataMapping','scaled');


% imshow(baseImg,[]);
contrastWindow = imcontrast(handles.axes1);
waitfor(contrastWindow);
newContrast = getDisplayRange(getimagemodel(findobj('parent',handles.axes1,'type','image')));
% baseImg(baseImg<newContrast(1)) = newContrast(1);
% baseImg(baseImg>newContrast(2)) = newContrast(2);
% baseImg = (baseImg-newContrast(1))./diff(newContrast);
% setappdata(handles.figure1,'displayImg',baseImg);
setappdata(handles.figure1,'newContrast',newContrast);

% currentColorMask = double(repmat(baseImg,[1,1,3]));
% currentColorMask = currentColorMask*0.8+coloredLabels.*0.2.*(1/255);
% set(imageHandle,'cdataMapping','direct');
% set(imageHandle,'cdata',currentColorMask);


% --- Executes on button press in showAll.
function showAll_Callback(hObject, eventdata, handles)
% hObject    handle to showAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
trackData=getappdata(handles.figure1,'trackOutput');
smoothWindow=str2double(get(handles.smoothingWindow,'String'));
startTime=str2double(get(handles.startTime,'string'));
normalizeFlag=get(handles.normalizeButton,'value');
timeStep=str2double(get(handles.timeStep,'string'));
minFrame=str2double(get(handles.minTime,'string'));
maxFrame=str2double(get(handles.maxTime,'string'));
minFrame=round(minFrame/timeStep);
maxFrame=round(maxFrame/timeStep);

trackData(trackData(:,7)<minFrame|trackData(:,7)>maxFrame,:)=[];
switch get(handles.plotChannel,'value')
    case 1
        output=trackData(:,4);
    case 2
        output=trackData(:,3);
        
    case 3
        output=trackData(:,4)./trackData(:,5);
end

output(trackData(:,end-1)<startTime)=nan;
output=normalizeRange(output);
%output=output/median(output);
nTracks=max(trackData(:,end));
nTime=max(trackData(:,end-1));
activityMat=zeros(nTracks,nTime);
for i=1:nTracks
    t=trackData((trackData(:,end)==i),end-1);
    a=output((trackData(:,end)==i));
    
    if normalizeFlag
a=(smooth(a/nanmedian(a),smoothWindow))/5-1;
a=a-min(a);
    end
    
    
    if  get(handles.interpMissing,'Value')
    if nnz(a~=0)>10
        [t,ib]=unique(t);
        a=a(ib);
a=interp1(t,a,min(t):max(t),'linear');
t=min(t):max(t);
    end
    end
        
    a=(smooth(a,smoothWindow));

    activityMat(i,t)=a;
    ngood(i)=nnz(trackData(:,end)==i);
end
activityMat(ngood/size(activityMat,2)<.5,:)=[];
%activityMat(activityMat==0)=nan;


setappdata(handles.figure1,'activityMat',activityMat);
cla(handles.axes2);
plotter(handles.slider1,eventdata);

% Hint: get(hObject,'Value') returns toggle state of showAll


% --- Executes on button press in alignmentSelect.
function alignmentSelect_Callback(hObject, eventdata, handles)
% hObject    handle to alignmentSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



[rpath,parent]=uigetfile('Y:\CommunalCode\3dbrain\');
registration=load([parent filesep rpath]);


setappdata(0,'registration',registration);


% --- Executes on slider movement.
function zSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function zSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in centerlineSelect.
function centerlineSelect_Callback(hObject, eventdata, handles)
% hObject    handle to centerlineSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%  load syncing data
dataFolder=uipickfiles('filterspec','E:');
dataFolder=dataFolder{1};
[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder);

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
hiResData=highResTimeTraceAnalysisTriangle3(dataFolder,imSize(1),imSize(2));
end

hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
[idxOut,bf2Hi]=flashTimeAlign2(bfFlashTime,hiResFlashTime);
f_hiResTime=fit(hiResFlashTime,bfFlashTime(bf2Hi),'poly1');
hiResData.frameTime=f_hiResTime(hiResData.frameTime);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

[idxOut,bf2fluor]=flashTimeAlign2(bfFlashTime,fluorFlashTime);
f_fluorTime=fit(fluorFlashTime,bfFlashTime(bf2fluor),'poly1');
fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

%works amazingly well, flash frames off by less than .05 seconds. 
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');


centerLineFile=dir([dataFolder filesep '*centerline*']);
centerLineFile={centerLineFile.name}';
if length(centerLineFile)>1
    centerlineFile=uipickfiles('FilterSpec',dataFolder);
    centerline=load(centerlineFile{1},'centerline');
    centerline=centerline.centerline;
else
centerline=load([dataFolder filesep centerLineFile{1}],'centerline');
centerline=centerline.centerline;
end
setappdata(handles.figure1,'bfIdxLookup',bfIdxLookup);
setappdata(handles.figure1,'centerline',centerline);



function maxTime_Callback(hObject, eventdata, handles)
% hObject    handle to maxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxTime as text
%        str2double(get(hObject,'String')) returns contents of maxTime as a double
timeStep=str2double(get(handles.timeStep,'string'));
maxFrame=str2double(get(hObject,'String'))/timeStep;

set(handles.slider1,'Max', maxFrame)

% --- Executes during object creation, after setting all properties.
function maxTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minTime_Callback(hObject, eventdata, handles)
% hObject    handle to minTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minTime as text
%        str2double(get(hObject,'String')) returns contents of minTime as a double

timeStep=str2double(get(handles.timeStep,'string'));

minFrame=str2double(get(hObject,'String'))/timeStep;
set(handles.slider1,'Min',minFrame )


% --- Executes during object creation, after setting all properties.
function minTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interpMissing.
function interpMissing_Callback(hObject, eventdata, handles)
% hObject    handle to interpMissing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interpMissing
showAll_Callback(handles.showAll, eventdata, handles)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, evnt, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
    if strcmp(evnt.Key,'rightarrow')|| strcmp(evnt.Key,'d')
            goForward_Callback(handles.slider1,evnt,handles);
            
        %Backward
        elseif strcmp(evnt.Key,'backspace') || strcmp(evnt.Key,'leftarrow')|| strcmp(evnt.Key,'a')
            goBack_Callback(handles.slider1,evnt,handles);
        %Up
        elseif  strcmp(evnt.Key,'uparrow')|| strcmp(evnt.Key,'w')
            set(handles.zSlider,'value',get(handles.zSlider,'Value')+1);
            plotter(handles.slider1,evnt);

        %Down
        elseif strcmp(evnt.Key,'downarrow')|| strcmp(evnt.Key,'s')
            set(handles.zSlider,'value',get(handles.zSlider,'Value')-1);
   plotter(handles.slider1,evnt);

    elseif strcmp(evnt.Key,'space')
            cursorNeuronSelect(handles.slider1,evnt)
        elseif strcmp(evnt.Key,'e');
            goForward_Callback(handles.slider1,evnt,handles);
            autoSelect(handles.slider1,evnt)
        elseif strcmp(evnt.Key,'1') 
            selectNeuron1_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt)

        elseif strcmp(evnt.Key,'2') 
            selectNeuron2_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'3')
            selectNeuron3_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt)
        elseif strcmp(evnt.Key,'4') 
            selectNeuron4_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'5') 
            selectNeuron5_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'6') 
            selectNeuron6_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'7')
            selectNeuron7_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'8') 
            selectNeuron8_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'9');
            selectNeuron9_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'0');
            selectNeuron10_Callback(handles.slider1,evnt,handles);
        elseif strcmp(evnt.Key,'-');
            selectNeuron11_Callback(handles.slider1,evnt,handles);
        elseif strcmp(evnt.Key,'=');
            selectNeuron12_Callback(handles.slider1,evnt,handles);
        elseif strcmp(evnt.Key,'z');
            selectNeuronN_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'x')
           selectNeuronM_Callback(handles.slider1,evnt,handles);
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'return');
            cursorNeuronSelect(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'shift');
            current=get(handles.channelSelect,'Value');
            set(handles.channelSelect,'Value',3-current);
             channelSelect_Callback(handles.channelSelect, evnt, handles)
        elseif strcmp(evnt.Key,'q');
            reclick(handles.slider1,evnt);
        elseif strcmp(evnt.Key,'c');
            exactNeuronSelect(handles.slider1,evnt)
        elseif strcmp(evnt.Key,'h')
            switchShow(handles,evnt)
        end

