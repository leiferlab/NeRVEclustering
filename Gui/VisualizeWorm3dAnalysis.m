function varargout = VisualizeWorm3dAnalysis(varargin)
% VISUALIZEWORM3DANALYSIS MATLAB code for VisualizeWorm3dAnalysis.fig
%      VISUALIZEWORM3DANALYSIS, by itself, creates a new VISUALIZEWORM3DANALYSIS or raises the existing
%      singleton*.
%
%      H = VISUALIZEWORM3DANALYSIS returns the handle to a new VISUALIZEWORM3DANALYSIS or the handle to
%      the existing singleton*.
%
%      VISUALIZEWORM3DANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZEWORM3DANALYSIS.M with the given input arguments.
%
%      VISUALIZEWORM3DANALYSIS('Property','Value',...) creates a new VISUALIZEWORM3DANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisualizeWorm3dAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisualizeWorm3dAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisualizeWorm3dAnalysis

% Last Modified by GUIDE v2.5 26-Jun-2017 12:04:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VisualizeWorm3dAnalysis_OpeningFcn, ...
    'gui_OutputFcn',  @VisualizeWorm3dAnalysis_OutputFcn, ...
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


% --- Executes just before VisualizeWorm3dAnalysis is made visible.
function VisualizeWorm3dAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualizeWorm3dAnalysis (see VARARGIN)

% Choose default command line output for VisualizeWorm3dAnalysis
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @plotSlide);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'SliderStep',[1,1]);

hlistener2=addlistener(handles.slider2,'ContinuousValueChange',...
    @plotter);
setappdata(handles.slider2,'hlistener',hlistener2);
set(handles.slider2,'SliderStep',[1,1]);

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

% UIWAIT makes VisualizeWorm3dAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualizeWorm3dAnalysis_OutputFcn(hObject, eventdata, handles)
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


%select Folder
display('Select an image from the stack');
mostRecent=getappdata(0,'mostRecent');
dataFolder=uipickfiles('filterspec',mostRecent,...
    'Prompt', 'Select Data Folder');
dataFolder=dataFolder{1};
setappdata(0,'mostRecent',fileparts(dataFolder));

%get CLstraight folder with the individual volume data and images
psFolder=dir([dataFolder filesep 'CLstraight*']);
imFolder=[dataFolder filesep psFolder(end).name];
imFiles=dir([imFolder filesep  '*.tif']);
isDigits=isstrprop(imFiles(1).name,'digit');
nDigits=sum(isDigits);
fileNameRoot=imFiles(1).name(1:find(isDigits,1,'first')-1);


display('Loading mat file with neuron coordinates');
pointStatsFile=[dataFolder filesep 'pointStatsNew.mat'];
dataMat=load(pointStatsFile);
fieldName=fieldnames(dataMat);
dataMat=getfield(dataMat,fieldName{1});
%save the pointStats file
setappdata(handles.figure1,'TrackData',dataMat);

display('Loading mat file with centerlines');
[centerline, offset]=loadCLBehavior(dataFolder);
CLdata.centerline=centerline;
CLdata.offset=offset;
%save centerline
setappdata(handles.figure1,'CLdata',CLdata);


display('Loading mat file with neural signals');
heatFile=[dataFolder filesep 'heatData.mat'];
heatData=load(heatFile);
%save centerline
setappdata(handles.figure1,'heatData',heatData);

%display the ethogram

% define ethogramColormap
fcolor=[0 1 0];%[27 158 119]/256;%green
bcolor=[1 0 0];%[217 95 2]/256;%red
turncolor=[0 0 1];%[117 112 179]/256;%blue
pausecolor=[255 217 50]/256;%yellow
ethocolormap=[bcolor;pausecolor;fcolor;turncolor];

try
    %this is for old version, will remove 
imagesc(heatData.ethoTrack','parent',handles.axes5)
catch me
    imagesc(heatData.behavior.ethogram','parent',handles.axes5)
end


caxis(handles.axes5,[-1 2]);
colormap(handles.axes5,ethocolormap);
axis(handles.axes5,'off');
hold(handles.axes5,'on');

%draw verticle line for current time
ethoPlot=plot(handles.axes5,[1 1], [ .5 1.5],'black');
setappdata(handles.figure1,'ethoPlot',ethoPlot);
hold(handles.axes5,'off');

imageInfo=imfinfo([imFolder filesep imFiles(1).name]);
setappdata(handles.figure1,'imFiles',imFiles);
setappdata(handles.figure1,'imFolder',imFolder);
setappdata(handles.figure1,'ndigits',nDigits);
setappdata(handles.figure1,'fileNameRoot',fileNameRoot);
%setting slider parameters
set(handles.slider1,'Min',1)

    set(handles.slider1,'Max',length(imFiles));

set(handles.slider1,'Value',1)

setappdata(handles.figure1,'currentFrame',1);
set(handles.slider2,'min',1);
set(handles.slider2,'max',length(imageInfo));
set(handles.slider2,'value',1);

[bfAll,~,hiResData]=tripleFlashAlign(dataFolder);
n_CL=size(centerline,3);
clTime=bfAll.frameTime(1:n_CL);
lookup=interp1(clTime,1:length(bfAll.frameTime),hiResData.frameTime(diff(hiResData.stackIdx)==1),'nearest',0);
size(lookup)
setappdata(handles.figure1,'lookup',lookup);
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


% --- Executes on selection change in channelSelect.
function channelSelect_Callback(hObject, eventdata, handles)
% hObject    handle to channelSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,eventdata);
% Hints: contents = cellstr(get(hObject,'String')) returns channelSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channelSelect


%pull the slice from the stack specified by the sliders
function baseImg= getImage(handles)

fileNameRoot=getappdata(handles.figure1,'fileNameRoot');
ndigits=getappdata(handles.figure1,'ndigits');
imFolder=getappdata(handles.figure1,'imFolder');

iImage=round(get(handles.slider1,'Value'));
iSlice=round(get(handles.slider2,'Value'));
imageName=[fileNameRoot num2str(iImage,['%6.' num2str(ndigits) 'd'])];


if get(handles.channelSelect,'Value')==2 || get(handles.channelSelect,'Value')==3
    dataFrame=getappdata(handles.figure1,'dataFrame');
    if any(dataFrame==iImage)
        ps=getappdata(handles.figure1,'ps');
    else
        ps=load([imFolder filesep 'pointStats' imageName(6:end) '.mat']);
        ps=ps.pointStats;
        setappdata(handles.figure1,'ps',ps);
    end
    baseMask=ps.baseImg(:,:,iSlice);
    contour(handles.axes1,baseMask,'black','LineWidth',.4);
end

if get(handles.channelSelect,'Value')==3
    baseImg=bwlabeln(ps.baseImg,6);
    baseImg=baseImg/max(baseImg(:));
    baseImg=baseImg(:,:,iSlice);
    return
elseif get(handles.channelSelect,'Value')==2
    baseImg=ps.baseImg;
    baseImg=baseImg(:,:,iSlice);
    return
end


if ~exist([imFolder filesep imageName '.tif'],'file')
    %% add some sort of error
    return
end
%read actual image
baseImg=double(imread([imFolder filesep imageName],'tif',...
    'Index',iSlice));

% apply contrast, may change how I do this
setappdata(handles.figure1,'baseImg',baseImg);

newContrast=getappdata(handles.figure1,'newContrast');
if isempty(newContrast)
    newContrast=[min(baseImg(:)),max(baseImg(:))];
end
baseImg(baseImg<newContrast(1)) = newContrast(1);
baseImg(baseImg>newContrast(2)) = newContrast(2);
baseImg = (baseImg-newContrast(1))./diff(newContrast);




%%%% get the current coordinates for the volume in the straightened
%%%% coordinate system, also get the requested IDs
function[ XYZcoord,pointID]=getPoints(handles)
displayRange=str2double(get(handles.displayRange,'String'));
iSlice=round(get(handles.slider2,'Value'));
iImage=round(get(handles.slider1,'Value'));
TrackData=getappdata(handles.figure1,'TrackData');
% plot points if they're there
trackRange={TrackData.stackIdx}';
trackRange(cellfun(@(x) isempty(x),trackRange))={-1};
iTrack=find(cell2mat(trackRange)==iImage);
%show only points that have been assigned an ID
XYZcoord=TrackData(iTrack).straightPoints;
%display
set(handles.FrameIdx,'string',['t=' num2str(TrackData(iTrack).stackIdx,'%6.2f')]);

%plot the type of points selected
switch get(handles.pointShowType,'Value')
    case 1
        pointID=TrackData(iTrack).trackIdx;
    case 2
        pointID=1:length(TrackData(iTrack).straightPoints);
    case 3
        pointID=TrackData(iTrack).matchIdx;
    case 4
        pointID=TrackData(iTrack).trackIdx;
        if ~isempty(heatData) %relabel with clustered labels
            pointID(~isnan(pointID))=heatData.cgIdxRev(pointID(~isnan(pointID)));
        end
    otherwise
        pointID=TrackData(iTrack).trackIdx;
end
pointID=pointID(:);

if isempty(pointID)
display('No Points Found!')
return
end

%only leave points that are in or close to the current slice being
%displayed
if ~get(handles.showAllPoints,'Value')
    showPoints=abs(XYZcoord(:,3)-iSlice)<=displayRange & ~isnan(pointID(:,end));
else
    showPoints=abs(XYZcoord(:,3)-iSlice)<=displayRange ;
end

XYZcoord=XYZcoord(showPoints,[2 1 3:end]);
pointID=pointID(showPoints);


function plotTrace= plotSignal(handles)
heatData=getappdata(handles.figure1,'heatData');
iImage=round(get(handles.slider1,'Value'));

%get target neuron
regionSelect=str2double(get(handles.trackNeuron,'String'));
%%% WILL NEED TO ADD THIS BACK IN WHEN FOLLOWING CLUSTER IDX
%regionSelect=heatData.cgIdxRev(regionSelect);

%plot the selected heatmap data if present
switch get(handles.plotChannel,'Value')
    case 1
        plotTrace=heatData.Ratio2(regionSelect,:);
    case 2
        plotTrace=heatData.R2(regionSelect,:);
    case 3
        plotTrace=heatData.G2(regionSelect,:);
    case 4
        plotTrace=heatData.rRaw(regionSelect,:);
    case 5
        plotTrace=heatData.gRaw(regionSelect,:);
end

switch get(handles.plotChannel2,'Value')
    case 2
        plotTrace2=heatData.Ratio2(regionSelect,:);
    case 3
        plotTrace2=heatData.R2(regionSelect,:);
    case 4
        plotTrace2=heatData.G2(regionSelect,:);
    case 5
        plotTrace2=heatData.rRaw(regionSelect,:);
    case 6
        plotTrace2=heatData.gRaw(regionSelect,:);
    otherwise
        plotTrace2=[];
end

%plot the selected heatmap data if present
if isempty(plotTrace) || regionSelect<0
    plotTrace=[];
    return
end

colorOrder='rygb';
%color ball based on behavior, or if no ethotrack, just use default
if isfield(heatData,'ethoTrack')
    currentcolor=colorOrder(heatData.ethoTrack(iImage)+2);
else
    currentcolor='b';
end

%if the neurontracked was changed, replot the signal, otherwise, just move
%the dot
 oldNeuron=getappdata(handles.figure1,'regionSelect');
 if ~any(oldNeuron==regionSelect)
 plot(handles.axes3,plotTrace);
 hold(handles.axes3,'on')
scatter(handles.axes3,iImage,plotTrace(iImage),currentcolor,'filled');
 hold(handles.axes3,'off');
 else
     tracePoint=findobj(handles.axes3,'type','Scatter');
     tracePoint(1).XData=iImage;
     tracePoint(1).YData=plotTrace(iImage);
 end
 
 if ~isempty(plotTrace2) && ~any(oldNeuron==regionSelect)
     hold(handles.axes3,'on')
     plot(handles.axes3,plotTrace2);
     scatter(handles.axes3,iImage,plotTrace2(iImage),currentcolor,'filled');
     hold(handles.axes3,'off');
 elseif ~isempty(plotTrace2)
     tracePoint(2).XData=iImage;
     tracePoint(2).YData=plotTrace(iImage);
 else
     tracePoint(2).XData=[];
     tracePoint(2).YData=[];
 end
 
 
%slide the window around the current time
 xRange=iImage+[-200,200];
 xRange=xRange-min(xRange(1)-1,1);
 xlim(handles.axes3,xRange);
 ylim(handles.axes3,nanmean( plotTrace(xRange(1):xRange(2)))...
     +[-.5 1]*(nanstd( plotTrace(xRange(1):xRange(2)))*4));

 
 %%% Also plot centerlines
CLdata=getappdata(handles.figure1,'CLdata');
if ~isempty(CLdata)
    
lookup=getappdata(handles.figure1,'lookup');
CLidx=lookup(iImage); %need to correct here to go from volIndex to ClIndex.
currentCL=CLdata.centerline(:,:,round(CLidx));
currentCL=bsxfun(@minus,currentCL,mean(currentCL));


CL_handle=findobj(handles.axes4,'type','Line');
head_handle=findobj(handles.axes4,'type','Scatter');
if isempty(CL_handle) && isempty(head_handle)
plot(handles.axes4,currentCL(:,2),currentCL(:,1),'linewidth',4);
hold(handles.axes4,'on')
scatter(handles.axes4,currentCL(1,2),currentCL(1,1),['o' currentcolor],'filled')
hold(handles.axes4,'off');
axis(handles.axes4,[-400 400 -400 400],'off');
else
    CL_handle.XData=currentCL(:,2);
    CL_handle.YData=currentCL(:,1);
    head_handle.XData=currentCL(1,2);
    head_handle.YData=currentCL(1,1);
    head_handle.MarkerFaceColor=currentcolor;
end

end


function plotter(hObject,eventdata)
% function does the actual plotting of image and points

%get globals and gui settings
handles=guidata(get(hObject,'Parent'));
regionSelect=str2double(get(handles.trackNeuron,'String'));

%get slider values
iImage=round(get(handles.slider1,'Value'));
if iImage>handles.slider1.Max
    iImage=handles.slider1.Max;
    handles.slider1.Value=handles.slider1.Max;
end

iSlice=round(get(handles.slider2,'Value'));
if isempty(iSlice)
    iSlice=1;
    set(handles.slider2,'Value',1);
end
set(handles.sliceIdx,'string',['z=' num2str(iSlice,'%6.2f')]);

%update line position on ethogram
ethoPlot=findobj(handles.axes5,'type','Line');
ethoPlot.XData=[iImage iImage];

%get images and points, plot the signal
baseImg= getImage(handles);
plotSignal(handles);
[ pointsi,pointID]=getPoints(handles);

%display the image by updating axes handle
im_handle=findobj(handles.axes1,'type','Image');
if isempty(im_handle)
    imagesc(baseImg,'parent',handles.axes1);
    
else
    set(im_handle,'CData',baseImg);
end
caxis(handles.axes1,[0 1]);


%delete previous text and scatter objects, slower, but this time is small
%compared to time to load images. 
scatterObjs=findobj(handles.axes1,'type','Scatter');
textObjs=findobj(handles.axes1,'type','Text');
delete(scatterObjs);
delete(textObjs);
if ~isempty(pointID)
%get the numbers of the points being displayed, create a cell array of them
hold(handles.axes1,'on')
pointIDstr=cellstr(num2str(pointID));
pointIDstr=cellfun(@(x) strrep(x,'   ','-'),pointIDstr,'uniform',0);
pointIDstr=cellfun(@(x) strrep(x,'NaN',''),pointIDstr,'uniform',0);

%plot points and text
scatter(handles.axes1,pointsi(:,1),pointsi(:,2),'xr');
text(pointsi(:,1),pointsi(:,2),pointIDstr(:),...
    'color','w','parent',handles.axes1);
end
% show the tracked point in green.
if any(pointID==regionSelect)
    nSelectIdx=pointID(:)==regionSelect;
    text(pointsi(nSelectIdx,1),pointsi(nSelectIdx,2),pointIDstr(nSelectIdx),...
    'color','g','parent',handles.axes1);
end

hold(handles.axes1,'off')



function startTime_Callback(hObject, eventdata, handles)
% hObject    handle to startTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 plotter(handles.slider1,eventdata);
% Hints: get(hObject,'String') returns contents of startTime as text
%        str2double(get(hObject,'String')) returns contents of startTime as a double


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


function moveFrame(handles, step)
%moves from current frame to next frame
%get slider values for image and slice
iImage=round(get(handles.slider1,'Value'));
iSlice=round(get(handles.slider2,'Value'));

if isempty(iSlice)
    iSlice=1;
    set(handles.slider2,'Value',1);
end

%if we're tracking a neuron, try to go to its slice
TrackData=getappdata(handles.figure1,'TrackData');


%get the neuron being tracked
neuronId=round(str2double(get(handles.trackNeuron,'String')));
% get the data for the neurons in the frame about to be shown
TrackDatai=TrackData(iImage+step);
%get the type of index being used for tracking
switch get(handles.pointShowType,'Value')
    case 1
pointID=TrackDatai.trackIdx;
    case 2
pointID=TrackDatai.pointIdx;
    case 3
pointID=TrackDatai.matchIdx;
    case 4
pointID=TrackDatai.trackIdx;
heatData=getappdata(handles.figure1,'heatData');
if ~isempty(heatData)
    pointID(~isnan(pointID))=heatData.cgIdxRev(pointID(~isnan(pointID)));
end
    otherwise
pointID=TrackDatai.trackIdx;
end

% get the neuron whose index matches the tracked neuron
neuronIdx=(pointID==neuronId);
neuronIdx=find(neuronIdx,1,'first');

if any(neuronIdx)
    set(handles.trackNeuron,'backgroundColor',[1 1 1]);
    %if its there, set the slice to be shown to be the position in the next
    %volume
    iSliceTemp=round(TrackDatai.straightPoints(neuronIdx,3));
    if ~isnan(iSliceTemp)
        iSlice=iSliceTemp;
    end
    
else
    %if the neuron is missing, show a red background on the tracked neuron
    set(handles.trackNeuron,'backgroundColor',[1 0 0]);
end

%set new slider positions
set(handles.slider1,'value',iImage+step);
set(handles.slider2,'value',iSlice);

plotter(handles.slider1)



% --- Executes on button press in goBack.
function goBack_Callback(hObject, eventdata, handles)
% hObject    handle to goBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FPS=getappdata(handles.figure1,'FPS');
handles=guidata(get(hObject,'Parent'));

moveFrame(handles, -FPS)


% --- Executes on button press in goForward.
function goForward_Callback(hObject, eventdata, handles)
% hObject    handle to goForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(get(hObject,'Parent'));
FPS=getappdata(handles.figure1,'FPS');

moveFrame(handles, FPS)

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
plotter(handles.slider1,eventdata);



% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider




function plotSlide(hObject,eventdata)

handles=guidata(get(hObject,'Parent'));
moveFrame(handles,0)

% --- Executes on button press in goForward2.
function goForward2_Callback(hObject, eventdata, handles)
% hObject    handle to goForward2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.slider1,'value',get(handles.slider1,'value')+1);
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




function movieTitle_Callback(hObject, eventdata, handles)
% hObject    handle to movieTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movieTitle as text
%        str2double(get(hObject,'String')) returns contents of movieTitle as a double


% --- Executes on button press in makeMovie.
function makeMovie_Callback(hObject, eventdata, handles)
% hObject    handle to makeMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state=get(hObject,'value');
ax1movie=[get(handles.movieTitle,'String') '1'];
ax2movie=[get(handles.movieTitle,'String') '2'];
startTime=str2double(get(handles.startTime,'String'));
imFolder=getappdata(handles.figure1,'imFolder');

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


% --------------------------------------------------------------------
function adjustContrast_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to adjustContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% display the base image, calculate the display image;
baseImg = getappdata(handles.figure1,'baseImg');
imageHandle = findobj(handles.axes1,'type','image');
storedImage = get(imageHandle);
set(imageHandle,'cdata',baseImg,'cdataMapping','scaled');


% imshow(baseImg,[]);
contrastWindow = imcontrast(handles.axes1);
waitfor(contrastWindow);
newContrast = getDisplayRange(getimagemodel(findobj(handles.axes1,'type','image')));
baseImg(baseImg<newContrast(1)) = newContrast(1);
baseImg(baseImg>newContrast(2)) = newContrast(2);
baseImg = (baseImg-newContrast(1))./diff(newContrast);
setappdata(handles.figure1,'displayImg',baseImg);
setappdata(handles.figure1,'newContrast',newContrast);



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
    
elseif strcmp(evnt.Key,'shift')
    %toggle between mask and image
    currentChannel=handles.channelSelect.Value;
    if currentChannel<3
        handles.channelSelect.Value=3-currentChannel;
    end
    plotter(handles.slider1);
    
end



% --- Executes on button press in selectPoints.
function selectPoints_Callback(hObject, eventdata, handles)
% hObject    handle to selectPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imFolder=getappdata(handles.figure1,'imFolder');

display('Select mat file with points');
matFile=uipickfiles('Filterspec',fileparts(imFolder),...
    'Prompt','Select the PointStatsNew file');
dataMat=load(matFile{1});
fieldName=fieldnames(dataMat);
dataMat=getfield(dataMat,fieldName{1});
setappdata(handles.figure1,'TrackData',dataMat);
dataFolder=fileparts(matFile{1});
[centerline, offset]=loadCLBehavior(dataFolder);
CLdata.centerline=centerline;
CLdata.offset=offset;
setappdata(handles.figure1,'CLdata',CLdata);
 [fiducialPoints,z2ImageIdxOffset]=loadFiducialPoints(dataFolder);

[bfTime,hiResFrameTime,hasPoints,bfRange,hiResRange,hasPointsTime,lookup]=dataTimeAlignment(dataFolder,fiducialPoints);

setappdata(handles.figure1,'lookup',lookup);

function displayRange_Callback(hObject, eventdata, handles)
% hObject    handle to displayRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayRange as text
%        str2double(get(hObject,'String')) returns contents of displayRange as a double




function trackNeuron_Callback(hObject, eventdata, handles)
% hObject    handle to trackNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trackNeuron as text
%        str2double(get(hObject,'String')) returns contents of trackNeuron as a double
plotter(handles.slider1);




% --- Executes on selection change in pointShowType.
function pointShowType_Callback(hObject, eventdata, handles)
% hObject    handle to pointShowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pointShowType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pointShowType



% --- Executes on button press in loadHeat.
function loadHeat_Callback(hObject, eventdata, handles)
% hObject    handle to loadHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');
heatFile=uipickfiles('filterspec',mostRecent,...
    'Prompt','Select the heatData.mat file');
heatData=load(heatFile{1});
setappdata(handles.figure1,'heatData',heatData);

% ethogramColormap
fcolor=[0 1 0];%[27 158 119]/256;
bcolor=[1 0 0];%[217 95 2]/256;
turncolor=[0 0 1];%[117 112 179]/256;
pausecolor=[255 217 50]/256;
ethocolormap=[bcolor;pausecolor;fcolor;turncolor];


imagesc(heatData.ethoTrack','parent',handles.axes5)
caxis(handles.axes5,[-1 2]);
colormap(handles.axes5,ethocolormap);
axis(handles.axes5,'off');
hold(handles.axes5,'on');

ethoPlot=plot(handles.axes5,[1 1], [ .5 1.5],'black');
setappdata(handles.figure1,'ethoPlot',ethoPlot);
hold(handles.axes5,'off');


% --- Executes on button press in selectPSfile.
function selectPSfile_Callback(hObject, eventdata, handles)
% hObject    handle to selectPSfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mostRecent=getappdata(0,'mostRecent');
pointStatsFile=uipickfiles('filterspec',mostRecent,...
    'Prompt', 'Select Data Folder');
pointStatsFile=pointStatsFile{1};
dataMat=load(pointStatsFile);

display('Loading mat file with neuron coordinates');
fieldName=fieldnames(dataMat);
dataMat=getfield(dataMat,fieldName{1});
%save the pointStats file
setappdata(handles.figure1,'TrackData',dataMat);


% --- Executes on button press in instructions.
function instructions_Callback(hObject, eventdata, handles)
% hObject    handle to instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({...
    '1. Select Folder and choose a fully analyzed BrainScanner Folder',...
    '2. The Pointstats file is loaded by default, if you do not want to use the PointStatsNew file, select a different one',...
    '3. Browse the video to get an idea of the straightening and tracking quality.',...
    '',...
    'Controls',...
    'asdw or arrow keys: forward and back in time, up and down in Z',...
    'Display Range: how many slices above and below a neuron should the neuron be displayed in.'...
    'Track Neuron: Which neuron to follow through time and show the trace of',...
    },...
    'Instructions');
