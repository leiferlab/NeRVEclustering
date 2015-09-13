function varargout = scan3Dtiffseries_wormPoints(varargin)
% SCAN3DTIFFSERIES_WORMPOINTS MATLAB code for scan3Dtiffseries_wormPoints.fig
%      SCAN3DTIFFSERIES_WORMPOINTS, by itself, creates a new SCAN3DTIFFSERIES_WORMPOINTS or raises the existing
%      singleton*.
%
%      H = SCAN3DTIFFSERIES_WORMPOINTS returns the handle to a new SCAN3DTIFFSERIES_WORMPOINTS or the handle to
%      the existing singleton*.
%
%      SCAN3DTIFFSERIES_WORMPOINTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAN3DTIFFSERIES_WORMPOINTS.M with the given input arguments.
%
%      SCAN3DTIFFSERIES_WORMPOINTS('Property','Value',...) creates a new SCAN3DTIFFSERIES_WORMPOINTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scan3Dtiffseries_wormPoints_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scan3Dtiffseries_wormPoints_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scan3Dtiffseries_wormPoints

% Last Modified by GUIDE v2.5 12-Jun-2015 13:38:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @scan3Dtiffseries_wormPoints_OpeningFcn, ...
    'gui_OutputFcn',  @scan3Dtiffseries_wormPoints_OutputFcn, ...
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


% --- Executes just before scan3Dtiffseries_wormPoints is made visible.
function scan3Dtiffseries_wormPoints_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scan3Dtiffseries_wormPoints (see VARARGIN)

% Choose default command line output for scan3Dtiffseries_wormPoints
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @goForward_Callback);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'SliderStep',[1,1]);

hlistener2=addlistener(handles.slider2,'ContinuousValueChange',...
    @plotSlide);
setappdata(handles.slider2,'hlistener',hlistener2);
set(handles.slider2,'SliderStep',[1,1]);


% 
% [rpath,parent]=uigetfile('Y:\CommunalCode\3dbrain\','Select Registration File');
% registration=load([parent filesep rpath]);
% 

%setappdata(0,'registration',registration);

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

% UIWAIT makes scan3Dtiffseries_wormPoints wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scan3Dtiffseries_wormPoints_OutputFcn(hObject, eventdata, handles)
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

    
imFolder=getappdata(0,'imFolder');


%select Folder
display('Select an image from the stack');

rawImFolder=uipickfiles('filterspec',imFolder);
rawImFolder=rawImFolder{1};

if ~isdir(rawImFolder)
    [rawImFolder, fileNameRoot]=fileparts(rawImFolder);
    isDigits=isstrprop(fileNameRoot,'digit');
    nDigits=sum(isDigits);
    fileNameRoot=fileNameRoot(~isDigits);  
end


imFiles=dir([rawImFolder filesep fileNameRoot '0*.tif']);

imageInfo=imfinfo([rawImFolder filesep imFiles(1).name]);
setappdata(0,'imFiles',imFiles);
setappdata(0,'rawImFolder',rawImFolder);
setappdata(0,'imFolder',rawImFolder);
setappdata(0,'ndigits',nDigits);
setappdata(0,'fileNameRoot',fileNameRoot);
%setting slider parameters
set(handles.slider1,'Min',1)

    set(handles.slider1,'Max',length(imFiles));

set(handles.slider1,'Value',1)

setappdata(handles.figure1,'currentFrame',1);
set(handles.slider2,'min',1);
set(handles.slider2,'max',length(imageInfo));
set(handles.slider2,'value',1);

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
handles=guidata(get(hObject,'Parent'));
currentFrame=getappdata(handles.figure1,'currentFrame');
dataFrame=getappdata(handles.figure1,'dataFrame');
timeStep=str2double(get(handles.timeStep,'string'));
smoothWindow=str2double(get(handles.smoothingWindow,'String'));
startTime=str2double(get(handles.startTime,'String'));
normalizeFlag=get(handles.normalizeButton,'value');
displayRange=str2double(get(handles.displayRange,'String'));
ndigits=getappdata(0,'ndigits');
fileNameRoot=getappdata(0,'fileNameRoot');
imFolder=getappdata(0,'imFolder');
iImage=round(get(handles.slider1,'Value'));
iSlice=round(get(handles.slider2,'Value'));
rawImFolder=getappdata(0,'rawImFolder');

imFiles=getappdata(0,'imFiles');
if isempty(iSlice)
    iSlice=1;
    set(handles.slider2,'Value',1);
end

TrackData=getappdata(handles.figure1,'TrackData');
imageName=[fileNameRoot num2str(iImage,['%6.' num2str(ndigits) 'd'])];
nIdx=(cellfun(@(x) str2double(x), cellstr(imageName')));
nIdx=find(~isnan(nIdx) & ~imag(nIdx));
set(handles.currentFolder,'string', imageName);
imageNumber=str2double(imageName(nIdx));

% try
% neuronId=round(str2double(get(handles.trackNeuron,'String')));
% pointsi=TrackData{iImage};
% neuronInfo=pointsi(pointsi(:,4)==neuronId,:);
% iSlice=round(neuronInfo(3));
% catch
% end


% plot points if they're there
trackRange={TrackData.stackIdx}';
trackRange(cellfun(@(x) isempty(x),trackRange))={-1};
iTrack=find(cell2mat(trackRange)==imageNumber);
set(handles.FrameIdx,'string',['t=' num2str(iTrack,'%6.2f')]);
if exist([rawImFolder filesep imageName '.tif'],'file') && any(iTrack)

%iTrack=iImage;
shapeVector='xso^d<>';
        %show only points that have been assigned an ID
pointsi=TrackData(iTrack).straightPoints;
switch get(handles.pointShowType,'Value')
    case 1
pointID=TrackData(iTrack).trackIdx;
    case 2
pointID=1:length(TrackData(iTrack).straightPoints);
    case 3
pointID=TrackData(iTrack).matchIdx;
    case 4
      pointID=  [TrackData(iTrack).trackIdx TrackData(iTrack).pointIdx];
    otherwise
pointID=TrackData(iTrack).trackIdx;
end
pointID=pointID(:);
   nselect=str2double(get(handles.trackNeuron,'String'));
       regionSelect=pointID==round(nselect);
 
% if any(regionSelect)
%     iSlice=round(pointsi(regionSelect,3));
%     
% end
set(handles.sliceIdx,'string',['z=' num2str(iSlice,'%6.2f')]);


%show image 
stackSize=get(handles.slider2,'max');


if mod(iImage,2)==1 && get(handles.flipOdd,'Value')
    imSlice=stackSize-iSlice+1;
elseif mod(iImage,2)==0 && get(handles.flipEven,'Value')
        imSlice=stackSize-iSlice+1;
else
    imSlice=iSlice;
end
 
baseImg=double(imread([rawImFolder filesep imageName],'tif',...
    'Index',imSlice));


%clear current axes
%arrayfun(@(x) delete(x),get(handles.axes1,'children'))

setappdata(handles.figure1,'baseImg',baseImg);
newContrast=getappdata(handles.figure1,'newContrast');
if isempty(newContrast)
    newContrast=[min(baseImg(:)),max(baseImg(:))];
end
baseImg(baseImg<newContrast(1)) = newContrast(1);
baseImg(baseImg>newContrast(2)) = newContrast(2);
baseImg = (baseImg-newContrast(1))./diff(newContrast);
hold(handles.axes1,'on')
objCounter=1;
while objCounter<length(handles.axes1.Children)
    if ~isa(handles.axes1.Children(objCounter),'matlab.graphics.primitive.Image')
        delete(handles.axes1.Children(objCounter))
    else
        objCounter=objCounter+1;
    end
    
    
    
    
end

if isempty(getappdata(handles.figure1,'ax'))
ax1=imagesc(baseImg,'parent',handles.axes1);
setappdata(handles.figure1,'ax',ax1);
else
    ax1=(getappdata(handles.figure1,'ax'));
  %  ax1=imagesc(baseImg,'parent',handles.axes1);

    set(ax1,'CData',baseImg);
setappdata(handles.figure1,'ax',ax1);
end

caxis(handles.axes1,[0 1]);

if get(handles.channelSelect,'Value')==3 || get(handles.channelSelect,'Value')==4

if dataFrame==currentFrame || isempty(dataFrame);
    ps=getappdata(handles.figure1,'ps');
else
    ps=load([rawImFolder filesep 'pointStats' imageName(6:end) '.mat']);
    ps=ps.pointStats;
    setappdata(handles.figure1,'ps',ps);
end
    baseMask=ps.baseImg(:,:,imSlice);

    contour(handles.axes1,baseMask,'black','LineWidth',.4);
end
 if get(handles.channelSelect,'Value')==4
    baseImg=bwlabeln(ps.baseImg,6);
    baseImg=baseImg(:,:,imSlice);
  set(ax1,'CData',baseImg);

 end
    if ~get(handles.showAllPoints,'Value')
showPoints=abs(pointsi(:,3)-iSlice)<=displayRange & ~isnan(pointID(:,end));
    else
   showPoints=abs(pointsi(:,3)-iSlice)<=displayRange ;
    end

pointsi=pointsi(showPoints,[2 1 3:end]);
pointID=pointID(showPoints);
    if isfield('TrackData','regionLabel')
 pointsRegion=TrackData(iTrack).regionLabel;
pointsRegion=pointsRegion(showPoints);
    else
        pointsRegion=ones(1,length(pointID));
    end
    
pointIDstr=cellstr(num2str(pointID));
pointIDstr=cellfun(@(x) strrep(x,'   ','-'),pointIDstr,'uniform',0);
pointIDstr=cellfun(@(x) strrep(x,'NaN',''),pointIDstr,'uniform',0);
pointsIdx=0:5;
for i=1:length(pointsIdx);
    regionSelect=pointsRegion==pointsIdx(i);
    if any(regionSelect)
scatter(handles.axes1,pointsi(regionSelect,1),pointsi(regionSelect,2),shapeVector(1));
text(pointsi(regionSelect,1),pointsi(regionSelect,2),pointIDstr(regionSelect),...
    'color','w','parent',handles.axes1);
    end

if any(pointID(regionSelect)==nselect)
    nSelectIdx=pointID(regionSelect)==nselect;
    
    text(pointsi(nSelectIdx,1),pointsi(nSelectIdx,2),pointIDstr(nSelectIdx),...
    'color','g','parent',handles.axes1);
    
    
end

end
if isfield(TrackData(iImage),'errIdx')
errIdx=TrackData(iImage).errIdx;
errIdx=errIdx(showPoints);
text(pointsi(errIdx,1),pointsi(errIdx,2),pointID(errIdx),...
    'color','r','parent',handles.axes1);
   
end
   
% 
%    text(pointsi(:,1),pointsi(:,2),pointID,...
%     'color','r','parent',handles.axes1); 
hold(handles.axes1,'off')
end
setappdata(handles.figure1,'dataFrame',currentFrame)



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

iImage=round(get(handles.slider1,'Value'));
iSlice=round(get(handles.slider2,'Value'));
set(handles.FrameIdx,'string',['t=' num2str(iImage,'%6.2f')]);
imFiles=getappdata(0,'imFiles');
if isempty(iSlice)
    iSlice=1;
    set(handles.slider2,'Value',1);
end

TrackData=getappdata(handles.figure1,'TrackData');
trackRange={TrackData.stackIdx}';
trackRange(cellfun(@(x) isempty(x),trackRange))={-1};
iTrack=find(cell2mat(trackRange)==iImage);

TrackData=getappdata(handles.figure1,'TrackData');
try
neuronId=round(str2double(get(handles.trackNeuron,'String')));
TrackDatai=TrackData(iTrack-FPS);

switch get(handles.pointShowType,'Value')
    case 1
pointID=TrackDatai.trackIdx;
    case 2
pointID=TrackDatai.pointIdx;
    case 3
pointID=TrackDatai.matchIdx;
    case 4
      pointID=  [TrackDatai.trackIdx TrackDatai.pointIdx];
    otherwise
pointID=TrackDatai.trackIdx;
end



neuronIdx=(pointID==neuronId);
if any(neuronIdx)
        set(handles.trackNeuron,'backgroundColor',[1 1 1]);

iSliceTemp=round(TrackDatai.straightPoints(neuronIdx,3));
if ~isnan(iSliceTemp)
iSlice=iSliceTemp;
end

else
    set(handles.trackNeuron,'backgroundColor',[1 0 0]);
    
end
catch
end


set(handles.slider1,'value',iImage-FPS);
set(handles.slider2,'value',iSlice);
setappdata(handles.figure1,'currentFrame',get(handles.slider1,'value'));
plotter(handles.slider1,eventdata)


% --- Executes on button press in goForward.
function goForward_Callback(hObject, eventdata, handles)
% hObject    handle to goForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(get(hObject,'Parent'));

FPS=getappdata(handles.figure1,'FPS');

iImage=round(get(handles.slider1,'Value'));
iSlice=round(get(handles.slider2,'Value'));
set(handles.FrameIdx,'string',['t=' num2str(iImage,'%6.2f')]);
imFiles=getappdata(0,'imFiles');
if isempty(iSlice)
    iSlice=1;
    set(handles.slider2,'Value',1);
end

TrackData=getappdata(handles.figure1,'TrackData');
trackRange={TrackData.stackIdx}';
trackRange(cellfun(@(x) isempty(x),trackRange))={-1};
iTrack=find(cell2mat(trackRange)==iImage);


TrackData=getappdata(handles.figure1,'TrackData');
try
neuronId=round(str2double(get(handles.trackNeuron,'String')));
TrackDatai=TrackData(iTrack+FPS);

switch get(handles.pointShowType,'Value')
    case 1
pointID=TrackDatai.trackIdx;
    case 2
pointID=TrackDatai.pointIdx;
    case 3
pointID=TrackDatai.matchIdx;
    case 4
      pointID=  [TrackDatai.trackIdx TrackDatai.pointIdx];
    otherwise
pointID=TrackDatai.trackIdx;
end



neuronIdx=(pointID==neuronId);
if any(neuronIdx)
        set(handles.trackNeuron,'backgroundColor',[1 1 1]);
iSliceTemp=round(TrackDatai.straightPoints(neuronIdx,3));
if ~isnan(iSliceTemp)
iSlice=iSliceTemp;
end
else
    set(handles.trackNeuron,'backgroundColor',[1 0 0]);
    
end
catch
end

set(handles.slider1,'value',iImage+FPS);
set(handles.slider2,'value',iSlice);

setappdata(handles.figure1,'currentFrame',get(handles.slider1,'value'));
plotter(handles.slider1,eventdata)



% --- Executes on button press in goUp.
function goUp_Callback(hObject, eventdata, handles)
% hObject    handle to goDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider2,'value',get(handles.slider2,'value')+1);
setappdata(handles.figure1,'currentSlice',get(handles.slider2,'value'));
plotter(handles.slider2,eventdata)




% --- Executes on button press in goDown.
function goDown_Callback(hObject, eventdata, handles)
% hObject    handle to goDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.slider2,'value',get(handles.slider2,'value')-1);
setappdata(handles.figure1,'currentSlice',get(handles.slider2,'value'));
plotter(handles.slider2,eventdata)






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
baseImg(baseImg<newContrast(1)) = newContrast(1);
baseImg(baseImg>newContrast(2)) = newContrast(2);
baseImg = (baseImg-newContrast(1))./diff(newContrast);
setappdata(handles.figure1,'displayImg',baseImg);
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
trackData=getappdata(0,'trackOutput');
smoothWindow=str2double(get(handles.smoothingWindow,'String'));
startTime=str2double(get(handles.startTime,'string'));
normalizeFlag=get(handles.normalizeButton,'value');


switch get(handles.plotChannel,'value')
    case 1
        output=trackData(:,4);
    case 2
        output=trackData(:,3);
      
    case 3
        output=trackData(:,3)./trackData(:,4);
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
    a=a(t>startTime);
    t=t(t>startTime);        
    a=(smooth(a,smoothWindow));
    if normalizeFlag
        a=normalizeRange(a);
    end
    
    activityMat(i,t)=a;
    
end
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


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, evnt, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



        %Forward
        if strcmp(evnt.Key,'rightarrow')|| strcmp(evnt.Key,'d')
            goForward_Callback(handles.slider1,evnt,handles);
            
        %Backward
        elseif strcmp(evnt.Key,'backspace') || strcmp(evnt.Key,'leftarrow')|| strcmp(evnt.Key,'a')
            goBack_Callback(handles.slider1,evnt,handles);
        %Up
        elseif  strcmp(evnt.Key,'uparrow')|| strcmp(evnt.Key,'w')
            goUp_Callback(handles.slider2,evnt,handles);
        %Down
        elseif strcmp(evnt.Key,'downarrow')|| strcmp(evnt.Key,'s')
            goDown_Callback(handles.slider2,evnt,handles);
  
        end
        


% --- Executes on button press in selectPoints.
function selectPoints_Callback(hObject, eventdata, handles)
% hObject    handle to selectPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


imFolder=getappdata(0,'imFolder');

display('Select mat file with points');
matFile=uipickfiles('Filterspec',fileparts(imFolder));
dataMat=load(matFile{1});
fieldName=fieldnames(dataMat);
dataMat=getfield(dataMat,fieldName{1});
setappdata(handles.figure1,'TrackData',dataMat);


% --- Executes on button press in flipOdd.
function flipOdd_Callback(hObject, eventdata, handles)
% hObject    handle to flipOdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flipOdd
plotter(hObject,eventdata);


% --- Executes on button press in flipEven.
function flipEven_Callback(hObject, eventdata, handles)
% hObject    handle to flipEven (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flipEven
plotter(hObject,eventdata);


% --- Executes on button press in showAllPoints.
function showAllPoints_Callback(hObject, eventdata, handles)
% hObject    handle to showAllPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showAllPoints
plotter(hObject,eventdata);



function displayRange_Callback(hObject, eventdata, handles)
% hObject    handle to displayRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayRange as text
%        str2double(get(hObject,'String')) returns contents of displayRange as a double


% --- Executes during object creation, after setting all properties.
function displayRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trackNeuron_Callback(hObject, eventdata, handles)
% hObject    handle to trackNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trackNeuron as text
%        str2double(get(hObject,'String')) returns contents of trackNeuron as a double


% --- Executes during object creation, after setting all properties.
function trackNeuron_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pointShowType.
function pointShowType_Callback(hObject, eventdata, handles)
% hObject    handle to pointShowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pointShowType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pointShowType


% --- Executes during object creation, after setting all properties.
function pointShowType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pointShowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
