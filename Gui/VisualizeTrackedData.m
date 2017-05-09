function varargout = VisualizeTrackedData(varargin)
% VISUALIZETRACKEDDATA MATLAB code for VisualizeTrackedData.fig
%      VISUALIZETRACKEDDATA, by itself, creates a new VISUALIZETRACKEDDATA or raises the existing
%      singleton*.
%
%      H = VISUALIZETRACKEDDATA returns the handle to a new VISUALIZETRACKEDDATA or the handle to
%      the existing singleton*.
%
%      VISUALIZETRACKEDDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZETRACKEDDATA.neuron6 with the given input arguments.
%
%      VISUALIZETRACKEDDATA('Property','Value',...) creates a new VISUALIZETRACKEDDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisualizeTrackedData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisualizeTrackedData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisualizeTrackedData

% Last Modified by GUIDE v2.5 05-May-2017 16:19:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VisualizeTrackedData_OpeningFcn, ...
    'gui_OutputFcn',  @VisualizeTrackedData_OutputFcn, ...
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


% --- Executes just before VisualizeTrackedData is made visible.
function VisualizeTrackedData_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualizeTrackedData (see VARARGIN)

% Choose default command line output for VisualizeTrackedData

%set up slider
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @moveFrame);
hlistenerz=addlistener(handles.zSlider,'ContinuousValueChange',...
    @moveFrame);
setappdata(handles.zSlider,'hlistener',hlistenerz);
set(handles.zSlider,'SliderStep',[1,1]);

setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'SliderStep',[1,1]);
setappdata(handles.figure1,'points',0)
setappdata(handles.figure1,'show',1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VisualizeTrackedData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualizeTrackedData_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SelectFolder.
function SelectFolder_Callback(~, eventdata, handles)
% hObject    handle to SelectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%select folers with data,
mostRecent=getappdata(0,'mostRecent');
display([...
    'For example, use this on VNC into tigress data and select the folder: '...
    '/tigress/LEIFER/PanNeuronal/testing_sets/BrainScanner20161031_111303/'...
    ])
if isempty(mostRecent)
    dataFolder=uipickfiles('Prompt', 'Select the BrainScanner Folder' );
else
    dataFolder=uipickfiles('filterspec', mostRecent,...
        'Prompt', 'Select the data Folder');
end
dataFolder=dataFolder{1};
mostRecent=dataFolder;
mostRecent=fileparts(mostRecent);

setappdata(0,'mostRecent',mostRecent);
set(handles.currentFolder,'String',dataFolder);
%load registration file if needed for split

setappdata(handles.figure1,'tifFlag',0);

registration=load([dataFolder filesep 'alignments.mat']);
if isfield(registration.alignments,'background')
    background=registration.alignments.background;
else
    background=0;
end
registration=registration.alignments.S2AHiRes;
setappdata(handles.figure1,'registration',registration)
setappdata(handles.figure1,'background',background);


if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder);
end

setappdata(handles.figure1,'hiResData',hiResData');
setappdata(handles.figure1,'imFiles',dataFolder);

%setting slider parameters
set(handles.slider1,'Min',1)
set(handles.slider1,'Value',1)
setappdata(handles.figure1,'cursorTarget', 1);
maxFrame=max(hiResData.stackIdx);
minZ=min(hiResData.Z);
maxZ=max(hiResData.Z);
set(handles.maxTime,'String',num2str(maxFrame));
set(handles.minTime,'String',num2str(1));
set(handles.slider1,'max',maxFrame);
setappdata(handles.figure1,'currentFrame',1);
set(handles.slider1,'value',1);
set(handles.zSlider,'min',minZ);
set(handles.zSlider,'max', maxZ);
set(handles.zSlider,'value',(maxZ+minZ)/2);

%%% prepare image object
baseImg=getImage(handles,1);
delete(findobj(handles.axes1,'type','image'))
imagesc(handles.axes1,baseImg)
hold(handles.axes1,'on')



%%% prepare text object
delete(findobj(handles.axes1,'type','text'))
for i=100:-1:1 %text objects are added in a stack, so count backwards
text( -20,-20,... %hide them off axis
    cellstr(num2str(i)),...  %premake the numbers 
    'VerticalAlignment','bottom',...
    'HorizontalAlignment','right',...
    'color',[1 1 1],...
    'fontsize',10,...
    'parent',handles.axes1);
end


%%% prepare scatter objects

delete(findobj(handles.axes1,'type','scatter'))
scatter(handles.axes1,[],[],'xr')
scatter(handles.axes1,[],[],'o')
hold(handles.axes1,'off')

%%% load fiducial points
fiducialFolder=[dataFolder filesep 'BotfFiducialPoints'];

%get all mat files in the fiducials folder. We used to have multiple users,
%but now there should only ever be one. 
fiducialFile=dir([fiducialFolder filesep '*.mat']);
fiducialFile=[fiducialFolder filesep fiducialFile(1).name];
if exist(fiducialFile,'file')
fiducialPoints=load(fiducialFile);
timeOffset=load([fiducialFolder filesep 'timeOffset']);
timeOffset=timeOffset.timeOffset;
setappdata(handles.figure1,'timeOffset',timeOffset)
fiducialPoints=fiducialPoints.fiducialPoints;
setappdata(handles.figure1,'fiducialPoints',fiducialPoints);
else
    statusWarning(handles.neuronCoordStatus,...
        'No neuron locations found! Has the botCheckCompiler run?',1)
end

%%% load heatmap data
heatDataFile=[dataFolder filesep 'heatData.mat'];
if exist(heatDataFile,'file')
heatData=load(heatDataFile);
setappdata(handles.figure1,'heatData',heatData);
if isfield(heatData,'behavior')
 statusWarning(handles.signalStatus,...
     'Neural Data Loaded',0);
else
 statusWarning(handles.signalStatus,...
     ' Warning: Signal found, but behavior missing!',0.5)
end
else
    statusWarning(handles.signalStatus,...
        'No neural activity found! Has the fiducialCropper run?',1)
end



% --- Executes on slider movement.
function slider1_Callback(~, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%Use the slider values to get the current frame to display
function [hiResIdx,iVolume]=getHiResIdx(handles)

%get slider positions, for time and Z
iVolume=round(get(handles.slider1,'Value'));
zPos=(get(handles.zSlider,'Value'));
if isempty(zPos)
    set(handles.zSlider,'Value',0)
    zPos=0;
end
zPos=zPos(1);
if isempty(zPos)
    zPos=1;
    set(handles.zSlider,'Value',1);
end

%pull data
offset=getappdata(handles.figure1,'timeOffset');
%for selection of dat files

hiResData=getappdata(handles.figure1,'hiResData');
vol_frame_list=getappdata(handles.figure1,'vol_frame_list');
zVoltages=getappdata(handles.figure1,'zVoltages');

% if we've moved to a new volume, re load the necessary voltage-frame
% relations and save them. 
if iVolume~=getappdata(handles.figure1,'currentFrame') || isempty(zVoltages)
    vol_frame_list=find(hiResData.stackIdx==iVolume);
    vol_frame_list=...
        vol_frame_list(vol_frame_list>(-offset)...
        & vol_frame_list<length(hiResData.stackIdx));
    zVoltages=hiResData.Z(vol_frame_list+offset);
    [zVoltages,~,ia]=unique(zVoltages);
    vol_frame_list=vol_frame_list(ia);
    [~,ib]=sort(zVoltages,'ascend');
    zVoltages=zVoltages(ib);
    vol_frame_list=vol_frame_list(ib);
    setappdata(handles.figure1,'vol_frame_list',vol_frame_list);
    setappdata(handles.figure1,'zVoltages',zVoltages);
end

zSlice=interp1(zVoltages,1:length(zVoltages),zPos,'nearest','extrap');
zVoltageOut=zVoltages(zSlice);
set(handles.zSlider,'Value',zVoltageOut);

hiResIdx=vol_frame_list(zSlice)+offset;
setappdata(handles.figure1,'currentHiResIdx',hiResIdx);
setappdata(handles.figure1,'currentFrame',iVolume);
set(handles.FrameIdx,'string',[num2str(iVolume,'%6.2f'), ...
    '  ' num2str(zVoltageOut)]);
setappdata(handles.figure1,'hiResIdx',hiResIdx);


% find the image at hiResIdx, reads data from handles to determine whether
% to show read or green image. 
function baseImg=getImage(handles,hiResIdx)

imFiles=getappdata(handles.figure1,'imFiles');
R=getappdata(handles.figure1,'registration');
background=getappdata(handles.figure1,'background');

%get image file ID, if its an error, reload it
Fid=getappdata(handles.figure1,'fileID');

sCMOSfile=dir([imFiles filesep '*.dat']);
sCMOSfile=sCMOSfile.name;
[row,col]=getdatdimensions(sCMOSfile);

if isempty(Fid)
    Fid=fopen([imFiles filesep sCMOSfile ] );
    setappdata(handles.figure1,'fileID',Fid);
elseif Fid<=0
    Fid=fopen([imFiles filesep sCMOSfile] );
    setappdata(handles.figure1,'fileID',Fid);
end

%move the pointer to the image and read it in
status=fseek(Fid,2*hiResIdx*row*col,-1);
fullImage=fread(Fid,row*col,'uint16',0,'l');
fullImage=(reshape(fullImage,row,col));

fullImage=fullImage-background;
fullImage(fullImage<0)=0;

% green and red image rectangles
rect1=R.rect1;
rect2=R.rect2;
t_concord=R.t_concord;
Rsegment=R.Rsegment;

%1 for red image, 0 for green
if get(handles.channelSelect,'Value')==1
    baseImg=fullImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
else
    activity=fullImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
    baseImg=imwarp(activity,t_concord,'OutputView',Rsegment);
end


%plots current frame in image window and plot window
function plotter(hObject,~)
%recover main handle
handles=guidata(get(hObject,'Parent'));



%get the frame index and the volume number
[hiResIdx,iVolume]=getHiResIdx(handles);

% pull up the image to display
baseImg=getImage(handles,hiResIdx);

% plot the actual image
ax1=findobj(handles.axes1,'type','Image');
if isempty(ax1)
    hold(handles.axes1,'on')
    imagesc(baseImg,'Parent',handles.axes1);
    axis(handles.axes1,'equal');
    hold(handles.axes1,'off')

else
    ax1.CData=baseImg;
end

newContrast=getappdata(handles.figure1,'newContrast');
if ~isempty(newContrast)
    caxis(handles.axes1, newContrast);
end



%remove current scatter and text objects from axes

pointPlotter(handles)

tracePlotter(handles);

function pointPlotter(handles)

%get flagged volumes and neurons
flagged_neurons=getappdata(handles.figure1,'flagged_neurons');
flagged_volumes=getappdata(handles.figure1,'flagged_volumes');

text_handles=findobj(handles.axes1,'type','text');
[hiResIdx,iVolume]=getHiResIdx(handles);

%get neuron points
fiducialPoints=getappdata(handles.figure1,'fiducialPoints');
if ~isempty(fiducialPoints) && length(fiducialPoints)>=iVolume
    currentFiducials=fiducialPoints{iVolume};
    
    if ~isempty(currentFiducials)
    setappdata(handles.figure1,'fiducials',currentFiducials);
   statusWarning(handles.neuronCoordStatus,'Neurons Present',0)
    else
        statusWarning(handles.neuronCoordStatus,'Warning, no neuron locations found!',.5)
        return
    end

else
        statusWarning(handles.neuronCoordStatus,'Warning, no neuron locations found!',.5)

    return
end

    % will will show the neuron currently being tracked, along with other
    % neurons that are either in the same plane or in a plane close by. 

currentTarget=str2double(get(handles.trackedNeuron,'String'));
currentTarget=round(currentTarget);
    % closeSlices are within +/- 2 frames of the currently shown frame,
    % these neurons will be shown in black
closePoints=cellfun(@(x) abs(x-hiResIdx)<3,currentFiducials(:,4),'uniform',0);
inSlicePoints=cellfun(@(x) x==hiResIdx,currentFiducials(:,4),'uniform',0);


%move text around, if the points are not close, move them off the screen
%rather than destroying them
for iNeuron=1:100
    if closePoints{iNeuron}
        text_handles(iNeuron).Position(1)=currentFiducials{iNeuron,1};
        text_handles(iNeuron).Position(2)=currentFiducials{iNeuron,2};
        
        if any(iVolume==flagged_volumes) || any(flagged_neurons==iNeuron)
            text_handles(iNeuron).Color=[.94 0 0];
        elseif iNeuron==currentTarget && ~all(text_handles(iNeuron).Color==[0 1 0])
            text_handles(iNeuron).Color=[0 1 0];
        elseif ~all(text_handles(iNeuron).Color==1)
            text_handles(iNeuron).Color=[1 1 1];
        end
    elseif any(text_handles(iNeuron).Position>0)
        text_handles(iNeuron).Position=[-10,-10,0];
    end
    
end

closePointsIds=cellfun(@(x) any(x),closePoints);

inSlicePointsIds=cellfun(@(x) any(x),inSlicePoints);
closeXY=cell2mat(currentFiducials(closePointsIds,1:2));
inSliceXY=cell2mat(currentFiducials(inSlicePointsIds,1:2));
scat_handles=findobj(handles.axes1,'type','scatter');
if  any(closePointsIds)
scat_handles(1).XData=closeXY(:,1);
scat_handles(1).YData=closeXY(:,2);
end
if any(inSlicePointsIds)
scat_handles(2).XData=inSliceXY(:,1);
scat_handles(2).YData=inSliceXY(:,2);
end



%plot the trace of the activity
function tracePlotter(handles)
target=str2double(get(handles.trackedNeuron,'String'));
target=round(target);
heatData=getappdata(handles.figure1,'heatData');

%if heatdata is missing, skip all of this
if isempty(heatData) || target<0
    return
end

%plot the selected data trace if present
plotType=get(handles.plotDisplay,'Value');
switch plotType
    case 1
        plotTrace=heatData.Ratio2(target,:);
    case 2
        plotTrace=heatData.R2(target,:);
    case 3
        plotTrace=heatData.G2(target,:);
    case 4
        plotTrace=heatData.rRaw(target,:);
    case 5
        plotTrace=heatData.gRaw(target,:);
end

iVolume=round(get(handles.slider1,'Value'));

if length(plotTrace)<iVolume
    statusWarning(handles.signalStatus, 'No signal found',.5)
    return
end

colorOrder='rygb';
%draw dot with colour based on behavior, if behavior is missing, use blue
if isfield(heatData,'behavior')
    currentBehavior=heatData.behavior.ethogram(iVolume);
    
    if isnan(currentBehavior)
        currentcolor=[ 1 1 1];% make occasional nans white
    else
    currentcolor=colorOrder(currentBehavior+2);
    end
else
    currentcolor='black';
end
oldPlotState=getappdata(handles.figure1,'oldPlot');
newPlotState=[target plotType];

%if either the neuron being tracked or the type of plot has changed, you
%need to reload the data, otherwise, just move the dot and slide the window
if ~all(ismember(oldPlotState,newPlotState)) || isempty(oldPlotState)
    plot(handles.axes3,plotTrace);
    hold(handles.axes3,'on')
    tracePoint=scatter(handles.axes3,...
        iVolume,plotTrace(iVolume),...
    currentcolor,'filled');
    setappdata(handles.figure1,'tracePoint',tracePoint);
    hold(handles.axes3,'off');
    setappdata(handles.figure1,'oldPlot',newPlotState)
else
    tracePoint=getappdata(handles.figure1,'tracePoint');
    tracePoint.YData=plotTrace(iVolume);
    tracePoint.XData=iVolume;
end
    statusWarning(handles.signalStatus, 'Signal found',0)

xlim(handles.axes3,[max(iVolume-100,1),min(length(plotTrace),iVolume+100)]);




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
step=str2double(get(handles.timeStep,'String'));
moveFrame(hObject,eventdata,-step);

%moves from the current frame to current frame + step. The program will try
%its best to keep the same Z level or stay on the same neuron. 
function moveFrame(hObject,eventdata,step)
if nargin==2
    step=0;
end

handles=guidata(get(hObject,'Parent'));
hiResData=getappdata(handles.figure1,'hiResData');

set(handles.slider1,'value',get(handles.slider1,'value')+step);
currentFrame=round(get(handles.slider1,'value'));
fiducialsAll=getappdata(handles.figure1,'fiducialPoints');
currentTarget=str2double(get(handles.trackedNeuron,'String'));
if length(fiducialsAll)>currentFrame
if size(fiducialsAll{currentFrame},1)>=currentTarget && size(fiducialsAll{currentFrame},2)>1
    newIdx=fiducialsAll{currentFrame}{currentTarget,4};
    newZ=hiResData.Z(newIdx);
else
    newZ=get(handles.zSlider,'value');
end
else
    newZ=get(handles.zSlider,'value');
end
%newZ=get(handles.zSlider,'value');

set(handles.zSlider,'value',min(newZ,get(handles.zSlider,'max')));


plotter(handles.slider1,eventdata)




% --- Executes on button press in goForward.
function goForward_Callback(hObject, eventdata, handles)
% hObject    handle to goForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
step=str2double(get(handles.timeStep,'String'));
moveFrame(hObject,eventdata,step)


% --- Executes on selection change in plotChannel.
function plotChannel_Callback(~, ~, handles)
% hObject    handle to plotChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,'eventdata');
% Hints: contents = cellstr(get(hObject,'String')) returns plotChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotChannel


% --- Executes during object creation, after setting all properties.
function plotChannel_CreateFcn(hObject, ~, ~)
% hObject    handle to plotChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




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


% --- Executes on button press in goUp.
function goUp_Callback(hObject, eventdata, handles)
% hObject    handle to goUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)set(handles.slider1,'value',get(handles.slider1,'value')+1);
zVoltages=getappdata(handles.figure1,'zVoltages');
currentZ=get(handles.zSlider,'value');
currentZ=currentZ(1);
if currentZ<max(zVoltages)
    newZ=zVoltages((find(zVoltages>currentZ,1,'first')));
else newZ=currentZ;
end
set(handles.zSlider,'value',min(newZ,get(handles.zSlider,'max')));
plotter(hObject,eventdata)

% --- Executes on button press in goDown.
function goDown_Callback(hObject, eventdata, handles)
% hObject    handle to goDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zVoltages=getappdata(handles.figure1,'zVoltages');
currentZ=get(handles.zSlider,'value');
if currentZ>min(zVoltages)
    newZ=zVoltages((find(zVoltages<currentZ,1,'last')));
else newZ=currentZ;
end
set(handles.zSlider,'value',min(newZ,get(handles.zSlider,'max')));
plotter(hObject,eventdata)




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
imageHandle = findobj('Parent',handles.axes1,'type','image');
set(imageHandle,'CData',imageHandle.CData,'cdataMapping','scaled');


contrastWindow = imcontrast(handles.axes1);
waitfor(contrastWindow);
newContrast = getDisplayRange(getimagemodel(findobj('parent',handles.axes1,'type','image')));
baseImg(baseImg<newContrast(1)) = newContrast(1);
baseImg(baseImg>newContrast(2)) = newContrast(2);
baseImg = (baseImg-newContrast(1))./diff(newContrast);
setappdata(handles.figure1,'displayImg',baseImg);
setappdata(handles.figure1,'newContrast',newContrast);

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
set(handles.slider1,'Value',minFrame);


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



% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, evnt, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
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
    goUp_Callback(handles.zSlider,evnt,handles);
    %Down
elseif strcmp(evnt.Key,'downarrow')|| strcmp(evnt.Key,'s')
    goDown_Callback(handles.zSlider,evnt,handles);
elseif strcmp(evnt.Key,'f');
    nextMissing_Callback(handles.slider1,evnt,handles)
elseif strcmp(evnt.Key,'z');
    previousMissing_Callback(handles.slider1,evnt,handles)
elseif strcmp(evnt.Key,'shift');
    current=get(handles.channelSelect,'Value');
    set(handles.channelSelect,'Value',3-current);
    channelSelect_Callback(handles.channelSelect, evnt, handles)
elseif strcmp(evnt.Key,'h')
    switchShow(handles,evnt)
end

function switchShow(handles,eventdata)

show=getappdata(handles.figure1,'show');
if isempty(show)
    show=1;
end
show=show+1;
show=mod(show,5);
setappdata(handles.figure1,'show',show);
switch getappdata(handles.figure1,'show')
    case 0
        display('Show none');
    case 1
        display('Show all points')
    case 2
        display('Show all points full name')
    case 3
        display('Show all points no labels');
    case 4
        display('Show only from current User');
end


plotter(handles.slider1,eventdata)

function statusWarning(hObject,statusString,fail)
if fail==1
    hObject.BackgroundColor=[.9 0 0];
elseif fail==0
        hObject.BackgroundColor=[0 1 0];
else
    hObject.BackgroundColor=[.9 .9 0];
end

hObject.String=statusString;



% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.figure1,'SelectionType'),'alt')
    %cursorNeuronSelect(handles.slider1,eventdata);
end



function maxIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to maxIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxIntensity as text
%        str2double(get(hObject,'String')) returns contents of maxIntensity as a double


% --- Executes during object creation, after setting all properties.
function maxIntensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in channelSelect.
function channelSelect_Callback(hObject, eventdata, handles)
% hObject    handle to channelSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channelSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channelSelect

plotter(handles.slider1,eventdata)


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


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
if eventdata.VerticalScrollCount>2
    for i=1:floor(eventdata.VerticalScrollCount/2)
        goUp_Callback(handles.zSlider,eventdata,handles);
    end
elseif eventdata.VerticalScrollCount<-2
    for i=1:floor(abs(eventdata.VerticalScrollCount/2))
        goDown_Callback(handles.zSlider,eventdata,handles);
    end
end



function trackedNeuron_Callback(hObject, eventdata, handles)
% hObject    handle to trackedNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trackedNeuron as text
%        str2double(get(hObject,'String')) returns contents of trackedNeuron as a double


% --- Executes during object creation, after setting all properties.
function trackedNeuron_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackedNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotDisplay.
function plotDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to plotDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotDisplay


% --- Executes during object creation, after setting all properties.
function plotDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadHeatData.
function loadHeatData_Callback(hObject, eventdata, handles)
% hObject    handle to loadHeatData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');
heatDataFile=uipickfiles('FilterSpec',mostRecent,...
    'Prompt', 'Select the heatData.mat file' );
heatDataFile=heatDataFile{1};
heatData=load(heatDataFile);
setappdata(handles.figure1,'heatData',heatData);

% --- Executes on selection change in pointShowType.
function pointShowType_Callback(~, ~, ~)
% hObject    handle to pointShowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pointShowType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pointShowType


% --- Executes during object creation, after setting all properties.
function pointShowType_CreateFcn(hObject, ~, ~)
% hObject    handle to pointShowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(~, ~, ~)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in flagNeuron.
function flagNeuron_Callback(hObject, ~, handles)
% hObject    handle to flagNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentTarget=str2double(get(handles.trackedNeuron,'String'));
flagged_neurons=getappdata(handles.figure1,'flagged_volumes');
flagged_neurons=[flagged_neurons, currentTarget];
setappdata(handles.figure1,'flagged_neurons',flagged_neurons);


% --- Executes on button press in flagTime.
function flagTime_Callback(hObject, eventdata, handles)
% hObject    handle to flagTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,iVolume]=getHiResIdx(handles);
flagged_volumes=getappdata(handles.figure1,'flagged_volumes');
flagged_volumes=[flagged_volumes, iVolume];
setappdata(handles.figure1,'flagged_volumes',flagged_volumes);


% --- Executes on button press in saveHeatMap.
function saveHeatMap_Callback(hObject, eventdata, handles)
% hObject    handle to saveHeatMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.currentFolder,'String',dataFolder);
flagged_volumes=getappdata(handles.figure1,'flagged_volumes');
flagged_neurons=getappdata(handles.figure1,'flagged_neurons');
%%% load heatmap data
heatDataFile=[dataFolder filesep 'heatData'];
if exist(heatDataFile,'file')
    save(heatDataFile,'flagged_volumes','flagged_neurons','-append')
else
    statusWarning(handles.signalStatus, ...
        'No neural activity found! Has the fiducialCropper run?')
end
