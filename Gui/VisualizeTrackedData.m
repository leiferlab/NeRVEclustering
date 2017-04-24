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

% Last Modified by GUIDE v2.5 17-Apr-2017 10:34:16

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
function VisualizeTrackedData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualizeTrackedData (see VARARGIN)

% Choose default command line output for VisualizeTrackedData

%set up slider
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @plotter);
hlistenerz=addlistener(handles.zSlider,'ContinuousValueChange',...
    @plotter);
setappdata(handles.zSlider,'hlistener',hlistenerz);
set(handles.zSlider,'SliderStep',[1,1]);
handles.timer = timer(...
    'ExecutionMode', 'fixedRate', ...       % Run timer repeatedly
    'Period', 60, ...                        % Initial period is 10 sec.
    'TimerFcn', {@updateData,hObject}); % Specify callback function

setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'SliderStep',[1,1]);
setappdata(handles.figure1,'points',0)
setappdata(handles.figure1,'show',1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VisualizeTrackedData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualizeTrackedData_OutputFcn(hObject, eventdata, handles)
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
mostRecent=getappdata(0,'mostRecent');
display([...
    'For example, use this on VNC into tigress data and select the folder: '...
    '/tigress/LEIFER/PanNeuronal/testing_sets/BrainScanner20161031_111303/'...
    ])
if isempty(mostRecent)
    dataFolder=uipickfiles('Prompt', 'Select the data Folder' );
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
setappdata(0,'registration',registration)
setappdata(0,'background',background);



if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    %only for 1200 by 600 image for now
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
minFrame=1;
set(handles.maxTime,'String',num2str(maxFrame));
set(handles.minTime,'String',num2str(minFrame));
set(handles.slider1,'max',maxFrame);
setappdata(handles.figure1,'currentFrame',1);
set(handles.slider1,'value',1);
set(handles.zSlider,'min',minZ);
set(handles.zSlider,'max', maxZ);
set(handles.zSlider,'value',(maxZ+minZ)/2);

%%% load fiducial points

fiducialFile=[dataFolder filesep 'BotfFiducialPoints'];
%fiducialData=load(fiducialFile);
timeOffset=load([fiducialFile filesep 'timeOffset']);
timeOffset=timeOffset.timeOffset;

setappdata(handles.figure1,'timeOffset',timeOffset)

%setappdata(handles.figure1,'fiducialsData', fiducialData2)
userList=dir([fiducialFile filesep '*.mat']); %get all user files
userList={userList.name};
userList=userList(~strcmp(userList,'timeOffset.mat'));
userList=cellfun(@(x) x(1:end-4),userList,'uniform',0);
fiducialPoints=load([fiducialFile filesep userList{1}]);
fiducialPoints=fiducialPoints.fiducialPoints;
setappdata(handles.figure1,'fiducialPoints',fiducialPoints);
stop(handles.timer)

start(handles.timer);


%%% load heatmap data
heatDataFile=[dataFolder filesep 'heatData'];
heatData=load(heatDataFile);
setappdata(handles.figure1,'heatData',heatData);





% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%plotter(handles.slider1,eventdata);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%setappdata(handles.figure1,'currentFrame',get(handles.slider1,'value'));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function plotter(hObject,eventdata)

%plots current frame in image window and plot window

handles=guidata(get(hObject,'Parent'));
iImage=round(get(handles.slider1,'Value'));
zPos=(get(handles.zSlider,'Value'));
zPos=zPos(1);
if isempty(zPos)
    zPos=1;
    set(handles.zSlider,'Value',1);
end

% minFrame=str2double(get(handles.minTime,'string'));
% maxFrame=str2double(get(handles.maxTime,'string'));

offset=getappdata(handles.figure1,'timeOffset');

imFiles=getappdata(handles.figure1,'imFiles');
%for selection of dat files

hiResData=getappdata(handles.figure1,'hiResData');
R=getappdata(0,'registration');
[row,col]=size(R.initialIm);
%row=1024;
%col=512;
FrameIdx=getappdata(handles.figure1,'FrameIdx');
zVoltages=getappdata(handles.figure1,'zVoltages');

if iImage~=getappdata(handles.figure1,'currentFrame') | isempty(zVoltages)
    
    FrameIdx=find(hiResData.stackIdx==iImage);%+offset;
    FrameIdx=FrameIdx(FrameIdx>(-offset) & FrameIdx<length(hiResData.stackIdx));
    zVoltages=hiResData.Z(FrameIdx);
    [zVoltages,~,ia]=unique(zVoltages);
    FrameIdx=FrameIdx(ia);
    [~,ib]=sort(zVoltages,'ascend');
    zVoltages=zVoltages(ib);
    FrameIdx=FrameIdx(ib);
    setappdata(handles.figure1,'FrameIdx',FrameIdx);
    setappdata(handles.figure1,'zVoltages',zVoltages);
    
    
end


stdPlot=hiResData.imSTD(FrameIdx+offset);
stdPlot=smooth(stdPlot,5);
[~,stdPeak]=max(stdPlot);
%zVoltages=zVoltages-zVoltages(stdPeak)+.2;
zSlice=interp1(zVoltages,1:length(zVoltages),zPos,'nearest','extrap');
zVoltageOut=zVoltages(zSlice);
set(handles.zSlider,'Value',zVoltageOut);
%    hiResIdx=metaData.iFrame(zSlice);
hiResIdx=FrameIdx(zSlice)+offset;
setappdata(handles.figure1,'currentHiResIdx',hiResIdx);
if getappdata(handles.figure1,'tifFlag')
    temp=getappdata(handles.figure1,'imageStack');
    baseImg=temp(:,:,hiResIdx);
else
    
Fid=getappdata(handles.figure1,'fileID');
sCMOSfile=dir([imFiles filesep '*.dat']);
sCMOSfile=sCMOSfile.name;
if isempty(Fid)
    Fid=fopen([imFiles filesep sCMOSfile ] );
    setappdata(handles.figure1,'fileID',Fid);
elseif Fid<=0;
    Fid=fopen([imFiles filesep sCMOSfile] );
    setappdata(handles.figure1,'fileID',Fid);
end

status=fseek(Fid,2*hiResIdx*row*col,-1);
temp=fread(Fid,row*col,'uint16',0,'l');
temp=(reshape(temp,row,col));

background=getappdata(0,'background');
temp=temp-background;
temp(temp<0)=0;

% fclose(Fid)
%     temp=pixelIntensityCorrection(temp);
%crop left and right regions
rect1=R.rect1;
rect2=R.rect2;
t_concord=R.t_concord;
Rsegment=R.Rsegment;
padRegion=R.padRegion;
worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

if get(handles.channelSelect,'Value')==1
    baseImg=worm; %red
else
    activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
    baseImg=imwarp(activity,t_concord,'OutputView',Rsegment);
end
end
setappdata(handles.figure1,'baseImg',baseImg);
setappdata(handles.figure1,'hiResIdx',hiResIdx);

% baseImg=pedistalSubtract(baseImg);
setappdata(handles.figure1,'currentFrame',iImage);
setappdata(0,'baseImg',baseImg)
%     figure
%     imagesc(smooth2a(baseImg,20,20)>5);
timeStep=str2double(get(handles.timeStep,'String'));
set(handles.FrameIdx,'string',[num2str(iImage*timeStep,'%6.2f') 's' ...
    '  ' num2str(zVoltageOut)]);
%    set(handles.FrameIdx,'string',stackName);

newContrast=getappdata(handles.figure1,'newContrast');
if isempty(newContrast)
    newContrast=[min(baseImg(:)),max(baseImg(:))];
end
%         baseImg(baseImg<newContrast(1)) = newContrast(1);
%         baseImg(baseImg>newContrast(2)) = newContrast(2);
%         baseImg = (baseImg-newContrast(1))./diff(newContrast);
hold(handles.axes1,'off')
ax1=findobj(handles.axes1,'type','Image');
if isempty(ax1)
    ax1=imagesc(baseImg,'Parent',handles.axes1);
else
    ax1.CData=baseImg;
end

caxis(handles.axes1, [newContrast]);
hold(handles.axes1,'on')
axis(handles.axes1,'equal');


colorOrder=get(gca,'colorOrder');
colorOrder=[colorOrder;colorOrder;colorOrder];
fiducialPoints=getappdata(handles.figure1,'fiducialPoints');
        textColor=[1 1 1] ;%textColor=colorOrder(iUser,:);
        textColor2=[0 0 0];
        textColor2=textColor;
    if ~isempty(fiducialPoints);
        %currentFiducialsAll=currentFiducialsAll.fiducialPoints;
        currentFiducials=fiducialPoints{iImage};
        
        setappdata(handles.figure1,'fiducials',currentFiducials);

        
        plotIdx=find(cell2mat((cellfun(@(x) ~isempty(x),currentFiducials(:,1),'uniformoutput',0))));
        currentPoints=cell2mat(currentFiducials(:,1:4));
    else
        currentPoints=[];
    end
    
    if ~isempty(currentPoints)
        
        %inSlice=interp1(zVoltages,1:length(zVoltages),currentPoints(:,4),'nearest','extrap');
        closeSlice=abs(currentPoints(:,4)-hiResIdx)<3;
        perfectSlice=currentPoints(:,4)==hiResIdx;
        currentTarget=str2double(get(handles.trackedNeuron,'String'));
        trackedPoint=currentPoints(currentTarget,:);
        delete(findobj(handles.axes1,'type','scatter'))
        delete(findobj(handles.axes1,'type','text'))
        if ~isempty(currentPoints)
         if closeSlice(currentTarget)
            text( trackedPoint(1), trackedPoint(2),num2str(currentTarget),...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right',...
                'color',[0 1 0],...
                 'fontsize',11,'parent',handles.axes1);
             scatter(handles.axes1,...
                 trackedPoint(1),...
                 trackedPoint(2),'black');
             if perfectSlice(currentTarget)
             scatter(handles.axes1,...
                 trackedPoint(1),...
                 trackedPoint(2),'xr');
             end
             closeSlice(currentTarget)=0;
             perfectSlice(currentTarget)=0;
             
             
         end
         
         
         switch getappdata(handles.figure1,'show')
                case 1
                    scatter(handles.axes1,currentPoints(closeSlice,1),currentPoints(closeSlice,2),'black');
                    scatter(handles.axes1,currentPoints(perfectSlice,1),currentPoints(perfectSlice,2),'xr');
                    text( currentPoints(closeSlice,1), currentPoints(closeSlice,2),[cellstr(num2str(plotIdx(closeSlice)))],'VerticalAlignment'...
                        ,'bottom', 'HorizontalAlignment','right','color',textColor,...
                        'fontsize',10,'parent',handles.axes1);
                    
                case 2
                    scatter(handles.axes1,currentPoints(closeSlice,1),currentPoints(closeSlice,2),'black');
                    scatter(handles.axes1,currentPoints(perfectSlice,1),currentPoints(perfectSlice,2),'xr');
                    text( currentPoints(closeSlice,1), currentPoints(closeSlice,2),...
                        cellfun(@(x) [user x],cellstr(num2str(plotIdx(closeSlice))),'uniform',0)...
                        ,'VerticalAlignment','bottom', ...
                        'HorizontalAlignment','right',...
                        'color',textColor,...
                        'fontsize',10,'parent',handles.axes1);
                case 3
                    cursorTarget=getappdata(handles.figure1,'cursorTarget');
                    cursorTarget=plotIdx==cursorTarget;
                    cursorTarget=cursorTarget & perfectSlice;
                    scatter(handles.axes1,currentPoints(cursorTarget,1),currentPoints(cursorTarget,2),'xr');
                case 4
                    if iUser ==get(handles.usersDropdown,'Value')
                    scatter(handles.axes1,currentPoints(closeSlice,1),currentPoints(closeSlice,2),'black');
                    scatter(handles.axes1,currentPoints(perfectSlice,1),currentPoints(perfectSlice,2),'xr');
                    text( currentPoints(closeSlice,1), currentPoints(closeSlice,2),[cellstr(num2str(plotIdx(closeSlice)))],'VerticalAlignment'...
                        ,'bottom', 'HorizontalAlignment','right','color',textColor,...
                        'fontsize',10,'parent',handles.axes1);  
                    end
         end
        end
        
    
    %allPoints=cat(2,allPoints,currentPoints);
end
%setappdata(handles.figure1,'allPoints',allPoints)
hold(handles.axes1,'off')
tracePlotter(handles,currentTarget);
drawnow;

function tracePlotter(handles,regionSelect)

heatData=getappdata(handles.figure1,'heatData');

%plot the selected heatmap data if present
if ~isempty(heatData) && regionSelect>0
    plotType=get(handles.plotDisplay,'Value');
 switch plotType
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

 iImage=round(get(handles.slider1,'Value'));
colorOrder='rygb';
%draw dot with colour based on behavior
currentBehavior=heatData.behavior.ethogram(iImage);
currentcolor=colorOrder(currentBehavior+2);
oldPlotState=getappdata(handles.figure1,'oldPlot');
newPlotState=[regionSelect plotType];


 if ~all(ismember(oldPlotState,newPlotState)) || isempty(oldPlotState)
 plot(handles.axes3,plotTrace);
 hold(handles.axes3,'on')
 tracePoint=scatter(handles.axes3,iImage,plotTrace(iImage),currentcolor,'filled');
 setappdata(handles.figure1,'tracePoint',tracePoint);
 hold(handles.axes3,'off');
 setappdata(handles.figure1,'oldPlot',newPlotState)
 else
     tracePoint=getappdata(handles.figure1,'tracePoint');
     tracePoint.YData=plotTrace(iImage);
     tracePoint.XData=iImage;
 end
 
 xlim(handles.axes3,[max(iImage-100,1),min(length(plotTrace),iImage+100)]);
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
step=str2double(get(handles.timeStep,'String'));
    moveFrame(hObject,eventdata,handles,-step);


function moveFrame(hObject,eventdata,handles,move)
set(handles.slider1,'value',get(handles.slider1,'value')+move);
currentSlide=get(handles.zSlider,'value');
currentFrame=round(get(handles.slider1,'value'));
previousFiducials=[];
currentFiducials=[];
fiducialsAll=getappdata(handles.figure1,'fiducialPoints');
steps=getappdata(handles.figure1,'stepCounter');

currentTarget=str2double(get(handles.trackedNeuron,'String'));
    
        currentPlotIdx=find(any(cell2mat((cellfun(@(x) ~isempty(x),...
            fiducialsAll{currentFrame},'uniformoutput',0))),2));
        
        currentFiducials=cell2mat(fiducialsAll{currentFrame}(:,1:4));
        
        
             if size(fiducialsAll{currentFrame},1)>=currentTarget && size(fiducialsAll{currentFrame},2)>1
                newZ=fiducialsAll{currentFrame}{currentTarget,3};
            else
                newZ=[];
            end

    
    %   fiducialPoints=getappdata(handles.figure1,'fiducials');
    %    currentFiducials=fiducialPoints{currentFrame};
    
    set(handles.zSlider,'value',min(newZ,get(handles.zSlider,'max')));
    

plotter(handles.slider1,eventdata)




% --- Executes on button press in goForward.
function goForward_Callback(hObject, eventdata, handles)
% hObject    handle to goForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%FPS=getappdata(handles.figure1,'FPS');
step=str2double(get(handles.timeStep,'String'));
moveFrame(hObject,eventdata,handles,step)


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



function timeOffset_Callback(hObject, eventdata, handles)
% hObject    handle to timeOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeOffset as text
%        str2double(get(hObject,'String')) returns contents of timeOffset as a double


% --- Executes during object creation, after setting all properties.
function timeOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeOffset (see GCBO)
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


% --- Executes on button press in loadFiducials.
function loadFiducials_Callback(hObject, eventdata, handles)
% hObject    handle to loadFiducials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imFiles=getappdata(handles.figure1,'imFiles');
if ~isdir(imFiles)
    parent=fileparts(imFiles);
else
    parent=imFiles;
end
fiducialFile=uipickfiles('filterspec',parent,...
   'Prompt', 'Select botFiducials Folder' );

fiducialFile=fiducialFile{1};
%fiducialData=load(fiducialFile);
timeOffset=load([fiducialFile filesep 'timeOffset']);
timeOffset=timeOffset.timeOffset;

setappdata(handles.figure1,'timeOffset',timeOffset)

%setappdata(handles.figure1,'fiducialsData', fiducialData2)
set(handles.currentFiducialFile,'String',fiducialFile)
userList=dir([fiducialFile filesep '*.mat']); %get all user files
userList={userList.name};
userList=userList(~strcmp(userList,'timeOffset.mat'));
userList=cellfun(@(x) x(1:end-4),userList,'uniform',0);
fiducialPoints=load([fiducialFile filesep userList{1}]);
fiducialPoints=fiducialPoints.fiducialPoints;
setappdata(handles.figure1,'fiducialPoints',fiducialPoints);
stop(handles.timer)

start(handles.timer);


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


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
