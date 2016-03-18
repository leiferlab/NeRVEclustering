function varargout = Worm3DFiducialPickerNetwork_weights(varargin)
% WORM3DFIDUCIALPICKERNETWORK_WEIGHTS MATLAB code for Worm3DFiducialPickerNetwork_weights.fig
%      WORM3DFIDUCIALPICKERNETWORK_WEIGHTS, by itself, creates a new WORM3DFIDUCIALPICKERNETWORK_WEIGHTS or raises the existing
%      singleton*.
%
%      H = WORM3DFIDUCIALPICKERNETWORK_WEIGHTS returns the handle to a new WORM3DFIDUCIALPICKERNETWORK_WEIGHTS or the handle to
%      the existing singleton*.
%
%      WORM3DFIDUCIALPICKERNETWORK_WEIGHTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORM3DFIDUCIALPICKERNETWORK_WEIGHTS.neuron6 with the given input arguments.
%
%      WORM3DFIDUCIALPICKERNETWORK_WEIGHTS('Property','Value',...) creates a new WORM3DFIDUCIALPICKERNETWORK_WEIGHTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Worm3DFiducialPickerNetwork_weights_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Worm3DFiducialPickerNetwork_weights_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Worm3DFiducialPickerNetwork_weights

% Last Modified by GUIDE v2.5 07-Oct-2015 04:38:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Worm3DFiducialPickerNetwork_weights_OpeningFcn, ...
    'gui_OutputFcn',  @Worm3DFiducialPickerNetwork_weights_OutputFcn, ...
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


% --- Executes just before Worm3DFiducialPickerNetwork_weights is made visible.
function Worm3DFiducialPickerNetwork_weights_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Worm3DFiducialPickerNetwork_weights (see VARARGIN)

% Choose default command line output for Worm3DFiducialPickerNetwork_weights

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

% UIWAIT makes Worm3DFiducialPickerNetwork_weights wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Worm3DFiducialPickerNetwork_weights_OutputFcn(hObject, eventdata, handles)
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


%select image folder
display('Select image folder, or image tif');
mostRecent=getappdata(0,'mostRecent');
if isempty(mostRecent)
    imFiles=uipickfiles();
else
    imFiles=uipickfiles('filterspec', mostRecent);
end
mostRecent=imFiles{1};
if ~isdir(mostRecent)
    mostRecent=fileparts(mostRecent);
end
setappdata(0,'mostRecent',mostRecent);
set(handles.currentFolder,'String',mostRecent);
%load registration file if needed for split
if strfind(imFiles{1}, '.tif');
    imageStack=stackLoad(imFiles{1});
    nImages=size(imageStack,3);
    hiResData.Z=1:nImages;
    hiResData.stackIdx=ones(1,nImages);
    hiResData.frametime=1:nImages;
    hiResData.imSTD=std(reshape(imageStack,[],nImages),1);
    setappdata(handles.figure1,'imageStack',imageStack);
    setappdata(handles.figure1,'tifFlag',1);
else
        setappdata(handles.figure1,'tifFlag',0);

if exist([mostRecent filesep 'alignments.mat'])
registration=load([mostRecent filesep 'alignments.mat']);
registration=registration.alignments.S2AHiRes;
setappdata(0,'registration',registration)
else
[rpath,parent]=uigetfile('Y:\CommunalCode\3dbrain\registration','Select Registration File');
registration=load([parent filesep rpath]);
end

setappdata(0,'registration',registration);

if exist([imFiles{1} filesep 'hiResData.mat'],'file')
    hiResData=load([imFiles{1} filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    %only for 1200 by 600 image for now
    hiResData=highResTimeTraceAnalysisTriangle4(imFiles{1},1200,600);
end

end
setappdata(handles.figure1,'hiResData',hiResData');
setappdata(handles.figure1,'imFiles',imFiles);

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






%function initialize(handles)

% --- Executes on button press in initialize.
function initialize_Callback(hObject, eventdata, handles)
% hObject    handle to initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');
hiResData=getappdata(handles.figure1,'hiResData');
fiducialPoints=cell(6,5);
fiducialPoints=repmat({fiducialPoints},max(hiResData.stackIdx),1);
setappdata(handles.figure1,'fiducials',fiducialPoints);
clickPoints=0;
contents = cellstr(get(handles.usersDropdown,'String'));
folderParent=uigetdir(mostRecent);
fiduicialFileName=[folderParent filesep datestr(now,'yyyymmddTHHMMSS') 'Fiducials' filesep];
set(handles.currentFiducialFile,'String',fiduicialFileName);
mkdir(fiduicialFileName);

for i=1:length(contents);
    save([fiduicialFileName filesep contents{i}],'fiducialPoints','clickPoints','-v6');
end
timeOffset=str2double(get(handles.timeOffset,'String'));
save([fiduicialFileName filesep 'timeOffset'],'timeOffset','-v6');
try
start(handles.timer);
catch 
end

%plotter(handles.slider1,eventdata);







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

offset=str2double(get(handles.timeOffset,'string'));

imFiles=getappdata(handles.figure1,'imFiles');
%for selection of dat files

hiResData=getappdata(handles.figure1,'hiResData');
R=getappdata(0,'registration');
%[row,col]=size(R.initialIm);
row=1200;
col=600;
imFiles=imFiles{1};
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
if isempty(Fid)
    Fid=fopen([imFiles filesep 'sCMOS_Frames_U16_1024x1024.dat'] );
    setappdata(handles.figure1,'fileID',Fid);
elseif Fid<=0;
    Fid=fopen([imFiles filesep 'sCMOS_Frames_U16_1024x1024.dat'] );
    setappdata(handles.figure1,'fileID',Fid);
end

status=fseek(Fid,2*hiResIdx*row*col,-1);
temp=fread(Fid,row*col,'uint16',0,'l');
temp=(reshape(temp,row,col));
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

ax1=imagesc(baseImg,'parent',handles.axes1);
set(ax1,'ButtonDownFcn',...
    'Worm3DFiducialPickerNetwork(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
caxis(handles.axes1, [newContrast]);
hold(handles.axes1,'on')
axis(handles.axes1,'equal');


%  B=bwboundaries(wormMask(:,:,zSlice));
%
%     for i=1:length(B)
%         b=B{i};
%         plot(handles.axes1,b(:,2),b(:,1),'b')
%     end
plot(handles.axes4,stdPlot,zVoltages);
ylim(handles.axes4,[get(handles.zSlider,'Min'),get(handles.zSlider,'Max')]);



%    fiducialData=getappdata(handles.figure1,'fiducialdata');

%fiducialFile=get(handles.currentFiducialFile,'String');
%   fiducialData=load(fiducialFile);
%  if isfield(fiducialData,'timeOffset');
%fiducialData = rmfield(fiducialData,'timeOffset');
%   end
%  fiducialData=structfun(@(x) x(iImage),fiducialData);
contents = cellstr(get(handles.usersDropdown,'String'));
fiducialFileName=get(handles.currentFiducialFile,'String');
colorOrder=get(gca,'colorOrder');
colorOrder=[colorOrder;colorOrder;colorOrder];
fiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
for iUser=1:length(contents);
    user=contents{iUser};
    if iUser==get(handles.usersDropdown,'Value')
        textColor=[1 1 1] ;%textColor=colorOrder(iUser,:);
        textColor2=[0 0 0];
    else
        textColor=colorOrder(iUser,:);
        textColor2=textColor;
    end
    if ~isempty(fiducialsAllUsers);
        currentFiducialsAll=fiducialsAllUsers.(user);
        %currentFiducialsAll=currentFiducialsAll.fiducialPoints;
        currentFiducials=currentFiducialsAll{iImage};
        
        if iUser==get(handles.usersDropdown,'Value');
            set(handles.DisplayIdx,'data',currentFiducials);
            setappdata(handles.figure1,'fiducials',currentFiducialsAll);
            
        end
        
        plotIdx=find(cell2mat((cellfun(@(x) ~isempty(x),currentFiducials(:,1),'uniformoutput',0))));
        currentPoints=cell2mat(currentFiducials(:,1:4));
        circleScatter=getappdata(handles.figure1,'circleScatter');
        if iUser==1
            delete(findobj(handles.axes2,'type','scatter'))
            delete(findobj(handles.axes2,'type','text'))
            cla(handles.axes2,'reset')
        end
    else
        currentPoints=[];
    end
    
    if ~isempty(currentPoints)
        hold(handles.axes2,'on')
        if iUser==get(handles.usersDropdown,'Value');
        randPos=(mod(10*sin(plotIdx.^(1.4)+iUser),5)-2)/2;
        circleScatter=scatter(handles.axes2,randPos,currentPoints(:,3),'x','MarkerEdgeColor',colorOrder(iUser,:));
        ylim(handles.axes2,ylim(handles.axes4))
        
        setappdata(handles.figure1,'circleScatter',circleScatter);
        
        end
        
        %inSlice=interp1(zVoltages,1:length(zVoltages),currentPoints(:,4),'nearest','extrap');
        closeSlice=abs(currentPoints(:,4)-hiResIdx)<str2double(get(handles.zSearch,'String'));
        perfectSlice=currentPoints(:,4)==hiResIdx;
        currentTarget=getappdata(handles.figure1,'cursorTarget');
        trackedPoint=currentPoints(currentTarget,:);
        if ~isempty(currentPoints)
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
                        ,'VerticalAlignment'...
                        ,'bottom', 'HorizontalAlignment','right','color',textColor,...
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
            text( trackedPoint(1), trackedPoint(2),num2str(currentTarget),'VerticalAlignment'...
                ,'bottom', 'HorizontalAlignment','right','color',[0 1 0],...
                'fontsize',10,'parent',handles.axes1);
            
            
        end
         if iUser ==get(handles.usersDropdown,'Value')
        circleLabel=text( randPos, currentPoints(:,3),...
            cellstr(num2str(plotIdx))...
            ,'VerticalAlignment'...
            ,'bottom', 'HorizontalAlignment','right','color',textColor2,...
            'fontsize',10,'parent',handles.axes2);
        setappdata(handles.figure1,'circleLabel',circleLabel);
         end
    end
    
    %allPoints=cat(2,allPoints,currentPoints);
end
plotCircle(handles)

%setappdata(handles.figure1,'allPoints',allPoints)
hold(handles.axes2,'off');

hold(handles.axes1,'off')
drawnow;



%     scat3=scatter3(handles.axes1,centroids(:,2),centroids(:,1),centroids(:,3),[],c);



%display circle on axis 2
function plotCircle(handles)
ylim(handles.axes2,[get(handles.zSlider,'Min'),get(handles.zSlider,'Max')]);

h=getappdata(handles.figure1,'circHandle');
if isempty(h)
    cla(handles.axes2)
    center=[0,.5];
    radius=1;
    h=viscircles(handles.axes2, center,radius);
    setappdata(handles.figure1,'circHandle',h);
    ylim(handles.axes2,[get(handles.zSlider,'Min'),get(handles.zSlider,'Max')]);
end


g=getappdata(handles.figure1,'lineHandle');
if isempty(g)
    
    hold(handles.axes2,'on')
    g=plot(handles.axes2,[-1,1], repmat(get(handles.zSlider,'Value'),1,2));
    hold(handles.axes2,'off');
    setappdata(handles.figure1,'lineHandle',g);
else
    if ishandle(g)
        set(g,'Ydata', repmat(get(handles.zSlider,'Value'),1,2));
    else
        hold(handles.axes2,'on')
        g=plot(handles.axes2,[-1,1], repmat(get(handles.zSlider,'Value'),1,2));
        hold(handles.axes2,'off');
        setappdata(handles.figure1,'lineHandle',g);
    end
    
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

    moveFrame(hObject,eventdata,handles,-1);


function moveFrame(hObject,eventdata,handles,move)
set(handles.slider1,'value',get(handles.slider1,'value')+move);
currentSlide=get(handles.zSlider,'value');
currentFrame=round(get(handles.slider1,'value'));
contents = cellstr(get(handles.usersDropdown,'String'));
previousFiducials=[];
currentFiducials=[];
tempFiducials=getappdata(handles.figure1,'fiducialsAllUsers');
steps=getappdata(handles.figure1,'stepCounter');
if  currentFrame>0
    
    for iUser=1:length(contents);
        user=contents{iUser};
        fiducialsAll=tempFiducials.(user);
        oldPlotIdx=find(any(cell2mat((cellfun(@(x) ~isempty(x),fiducialsAll{currentFrame-move},'uniformoutput',0))),2));
        currentPlotIdx=find(any(cell2mat((cellfun(@(x) ~isempty(x),fiducialsAll{currentFrame},'uniformoutput',0))),2));
        
        overlap=intersect(currentPlotIdx,oldPlotIdx);
        
        currentFiducials=cat(1,currentFiducials,fiducialsAll{currentFrame}(overlap,1:4));
        
        
        previousFiducials=cat(1,previousFiducials,fiducialsAll{currentFrame-move}(overlap,1:4));
        
        
        if iUser==get(handles.usersDropdown,'Value')
            currentTarget=getappdata(handles.figure1,'cursorTarget');
            % should always be true
            if size(fiducialsAll{currentFrame},1)>=currentTarget && size(fiducialsAll{currentFrame},2)>1
                newZ=fiducialsAll{currentFrame}{currentTarget,3};
            else
                newZ=[];
            end
        end
        
    end
    
    %   fiducialPoints=getappdata(handles.figure1,'fiducials');
    %    currentFiducials=fiducialPoints{currentFrame};
    
    if isempty(newZ)
        try
            if size(previousFiducials,1)>4
                oldzVoltages=(cell2mat(previousFiducials(:,3)));
                currentzVoltages=(cell2mat(currentFiducials(:,3)));
                p=polyfit(oldzVoltages,currentzVoltages,1);
                newZ=p(2)+p(1)*currentSlide;
            else
                newZ=currentSlide;
            end
        catch
            newZ=currentSlide;
        end
        
    end
    set(handles.zSlider,'value',min(newZ,get(handles.zSlider,'max')));
    
end
if steps>10
    steps=1;
else
    steps=steps+1;
end

plotter(handles.slider1,eventdata)




% --- Executes on button press in goForward.
function goForward_Callback(hObject, eventdata, handles)
% hObject    handle to goForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%FPS=getappdata(handles.figure1,'FPS');

moveFrame(hObject,eventdata,handles,1)


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
%cursorNeuronSelect(hObject,eventdata)
%displayIdx=get(handles.DisplayIdx,'data');
cursorValue=str2double(get(hObject,'String'));
setappdata(handles.figure1,'cursorTarget',cursorValue);
buttonClear(handles)
set(hObject,'BackgroundColor',[1 0 0 ]);

% --- Executes on button press in selectNeuron2.
function selectNeuron2_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursorValue=str2double(get(hObject,'String'));
setappdata(handles.figure1,'cursorTarget',cursorValue);
buttonClear(handles)
set(hObject,'BackgroundColor',[1 0 0 ]);

% --- Executes on button press in selectNeuron3.
function selectNeuron3_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursorValue=str2double(get(hObject,'String'));
setappdata(handles.figure1,'cursorTarget',cursorValue);
buttonClear(handles)
set(hObject,'BackgroundColor',[1 0 0 ]);


% --- Executes on button press in selectNeuron4.
function selectNeuron4_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursorValue=str2double(get(hObject,'String'));
setappdata(handles.figure1,'cursorTarget',cursorValue);
buttonClear(handles)
set(hObject,'BackgroundColor',[1 0 0 ]);


% --- Executes on button press in selectNeuron5.
function selectNeuron5_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursorValue=str2double(get(hObject,'String'));
setappdata(handles.figure1,'cursorTarget',cursorValue);
buttonClear(handles)
set(hObject,'BackgroundColor',[1 0 0 ]);


% --- Executes on button press in selectNeuron6.
function selectNeuron6_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursorValue=str2double(get(hObject,'String'));
setappdata(handles.figure1,'cursorTarget',cursorValue);
buttonClear(handles)
set(hObject,'BackgroundColor',[1 0 0 ]);


% --- Executes on button press in selectNeuron7.
function selectNeuron7_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',7);
%cursorNeuronSelect(hObject,eventdata)

% --- Executes on button press in selectNeuron8.
function selectNeuron8_Callback(hObject, eventdata, handles)
% hObject    handle to selectNeuron8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1,'cursorTarget',8);
%cursorNeuronSelect(hObject,eventdata)


function cursorNeuronSelect(hObject,eventdata)

handles=guidata(get(hObject,'Parent'));


[xselect,yselect]=ginput(1);
windowSearch=str2double(get(handles.xySearch,'String'));
zSearch=str2double(get(handles.zSearch,'String'));
cornersX=xselect+[windowSearch windowSearch -windowSearch -windowSearch windowSearch];
cornersY=yselect+[windowSearch -windowSearch -windowSearch windowSearch windowSearch];


inputData(hObject,xselect,yselect,windowSearch,zSearch,2)
setappdata(handles.figure1,'lastClick',[xselect yselect]);
hold(handles.axes1,'on')
plot(handles.axes1,cornersX,cornersY,'g')
hold(handles.axes1,'off');
drawnow

function exactNeuronSelect(hObject,eventdata)

handles=guidata(get(hObject,'Parent'));
[xselect,yselect]=ginput(1);
windowSearch=1;
zSearch=0;
cornersX=xselect+[windowSearch windowSearch -windowSearch -windowSearch windowSearch];
cornersY=yselect+[windowSearch -windowSearch -windowSearch windowSearch windowSearch];


inputData(hObject,xselect,yselect,windowSearch,zSearch,2)
setappdata(handles.figure1,'lastClick',[xselect yselect]);
hold(handles.axes1,'on')
plot(handles.axes1,cornersX,cornersY,'g')
hold(handles.axes1,'off');



function autoSelect(hObject,eventdata)
handles=guidata(get(hObject,'Parent'));
iFrame=getappdata(handles.figure1,'currentFrame');
lastFrame=get(handles.slider1,'Max');
oldFrameIdx=str2double(get(handles.refIdx,'String'));
if isempty(getappdata(handles.figure1,'cruiseStartFrame'))
    oldFrameIdx=iFrame+((oldFrameIdx:-2:-oldFrameIdx));
else
    oldFrameIdx=getappdata(handles.figure1,'cruiseStartFrame')+(-1:-1:-oldFrameIdx);
end

oldFrameIdx=oldFrameIdx(oldFrameIdx>0 & oldFrameIdx<lastFrame);

if get(handles.compareOnlyPrevious,'Value')
    oldFrameIdx=oldFrameIdx(oldFrameIdx<iFrame);
    
end
 oldFrameIdx=oldFrameIdx(oldFrameIdx~=iFrame);
contents = cellstr(get(handles.usersDropdown,'String'));
fiducialFileName=get(handles.currentFiducialFile,'String');
plotIdx=getappdata(handles.figure1,'cursorTarget');
userIdx=get(handles.usersDropdown,'Value');
%find best match frame for all users and all points


tempFiducials=getappdata(handles.figure1,'fiducialsAllUsers');

%search for previous data with similar distance matrix
bestOverlap=0;bestDistance=Inf;
currentAll=structfun(@(x) x{iFrame},tempFiducials,'uniform',0);
currentAll=struct2cell(currentAll);
%cellfun inside a cell fun to find index of fiducials... trust me it
%works
currentFiducialsIdx=cellfun(@(Y) ((any(cell2mat(cellfun(@(x) ~isempty(x),Y,...
    'uniformoutput',0)),2))),currentAll,'uniformOutput',0)  ;

currentFiducialsIdx{userIdx}(plotIdx)=0;
currentFiducialsLength=cellfun(@(x) length(x),currentFiducialsIdx);
for iCounter=1:length(oldFrameIdx)
    
    iRefFrame=oldFrameIdx(iCounter);
    
    oldFiducials=structfun(@(x) x{iRefFrame},tempFiducials,'uniform',0);
    oldFiducials=struct2cell(oldFiducials);
  
oldFiducialIdx=cellfun(@(Y) ((any((cellfun(@(x) ~isempty(x),Y(:,1))),2)))...
    ,oldFiducials,'uniformOutput',0)  ;
     oldFiducialsLength=cellfun(@(x) length(x),oldFiducialIdx);
    oldFiducialsLength=min(oldFiducialsLength,currentFiducialsLength);
    oldFiducialsLength=num2cell(oldFiducialsLength); 
    overlapIdx=cellfun(@(a,b,c) find(and(a(1:c),b(1:c))'),oldFiducialIdx,currentFiducialsIdx,oldFiducialsLength,'uniform',0);
    overlapN=sum(cellfun(@(x) length(x), overlapIdx));
    if any(cell2mat(oldFiducialIdx)) && length(oldFiducialIdx{userIdx})>=(plotIdx) && overlapN>2
        if oldFiducialIdx{userIdx}(plotIdx)
            controlPoint=cell2mat(oldFiducials{userIdx}(plotIdx,1:4));
        movingPointstemp=(cellfun(@(a,b) (a(b,1:5)'),oldFiducials,overlapIdx,'uniform',0));
        movingPointstemp=cell2mat([movingPointstemp{:}]');
        masterPointstemp=(cellfun(@(a,b) (a(b,1:5)'),currentAll,overlapIdx,'uniform',0));
        masterPointstemp=cell2mat([masterPointstemp{:}]');

        d2controlPoint=pdist2(controlPoint(:,1:3),movingPointstemp(:,1:3));
        d2controlPoint=sqrt(bsxfun(@plus,d2controlPoint.^2,d2controlPoint'.^2));

        d2controlPoint=d2controlPoint(tril(true(overlapN),-1))';
        
          d2Weights=1./d2controlPoint;
        movingDmat=pdist(movingPointstemp(:,1:3));
        masterDmat=pdist(masterPointstemp(:,1:3));
        overlap=sum(cell2mat(cellfun(@(x) nnz(x), overlapIdx,'uniformoutput',0)));
        dist=sum(d2Weights.*sqrt((movingDmat-masterDmat).^2));
        
        
        if overlap>=bestOverlap  && dist<=bestDistance && (oldFiducialIdx{userIdx}(plotIdx))
            bestRef=iRefFrame;
            bestOverlap=overlap;
            bestDistance=dist;
            movingPoints=movingPointstemp(:,[1 2 4]);
            masterPoints=masterPointstemp(:,[1 2 4]);
            plotIdxPoint=oldFiducials{userIdx}(plotIdx,:);
        end
        end
    end

end
dbquitdisplay(['Frame ' num2str(iFrame) ' matches ' num2str(bestRef)]);

plotIdxPoint=cell2mat(plotIdxPoint([1 2 4]));
movingFloor=min(movingPoints(:,3))-5;
masterFloor=min(masterPoints(:,3))-5;
movingPoints(:,end)=movingPoints(:,end)-movingFloor;
masterPoints(:,end)=masterPoints(:,end)-masterFloor;
plotIdxPoint(:,end)=plotIdxPoint(:,end)-movingFloor;
%
% [newEstimatePoint(:,1),newEstimatePoint(:,2),newEstimatePoint(:,3)]...
%     =transformPointsInverse(tform,plotIdxPoint(:,1),...
%     plotIdxPoint(:,2),plotIdxPoint(:,3));
% Fx=scatteredInterpolant(movingPoints(:,1),movingPoints(:,2),movingPoints(:,3),masterPoints(:,1));
% Fy=scatteredInterpolant(movingPoints(:,1),movingPoints(:,2),movingPoints(:,3),masterPoints(:,2));
% %Fz=scatteredInterpolant(movingPoints(:,1),movingPoints(:,2),movingPoints(:,3),masterPoints(:,3));
%
% Fx=scatteredInterpolant(movingPoints(:,1),movingPoints(:,2),masterPoints(:,1));
% Fy=scatteredInterpolant(movingPoints(:,1),movingPoints(:,2),masterPoints(:,2));
% %Fz=scatteredInterpolant(movingPoints(:,1),movingPoints(:,2),movingPoints(:,3));,masterPoints(:,3));
%
% newEstimatePoint=[Fx(plotIdxPoint(1),plotIdxPoint(2))...
%     Fy(plotIdxPoint(1),plotIdxPoint(2))];
tform = makeAffine3d(movingPoints, masterPoints);


[movingPoints(:,1),movingPoints(:,2),movingPoints(:,3)]...
    =transformPointsForward(tform,movingPoints(:,1),...
    movingPoints(:,2),movingPoints(:,3));
[plotIdxPoint(:,1),plotIdxPoint(:,2),plotIdxPoint(:,3)]...
    =transformPointsForward(tform,plotIdxPoint(:,1),...
    plotIdxPoint(:,2),plotIdxPoint(:,3));


newEstimatePoint=tpswarp3points(movingPoints,masterPoints,plotIdxPoint);
newEstimatePoint(3)=newEstimatePoint(3)+masterFloor;

%newZestimate=Fz(plotIdxPoint(1),plotIdxPoint(2))

% newEstimatePoint(1)=Fx(plotIdxPoint(1),plotIdxPoint(2),plotIdxPoint(3));
%newEstimatePoint(2)=Fy(plotIdxPoint(1),plotIdxPoint(2),plotIdxPoint(3));
%newEstim)atePoint(3)=Fz(plotIdxPoint(1),plotIdxPoint(2),plotIdxPoint(3));
hiResIdx=round(newEstimatePoint(3));
frameIdx=getappdata(handles.figure1,'FrameIdx');

hiResIdx=max(hiResIdx,min(frameIdx)+ str2double(get(handles.timeOffset,'String')));
hiResIdx=min(hiResIdx,max(frameIdx)+ str2double(get(handles.timeOffset,'String')));
setappdata(handles.figure1,'currentHiResIdx',hiResIdx);


windowSearch=str2double(get(handles.xySearch,'String'));
zSearch=str2double(get(handles.zSearch,'String'));

cornersX=newEstimatePoint(:,1)+[windowSearch windowSearch -windowSearch -windowSearch windowSearch];
cornersY=newEstimatePoint(:,2)+[windowSearch -windowSearch -windowSearch windowSearch windowSearch];

inputData(hObject,newEstimatePoint(:,1),newEstimatePoint(:,2),windowSearch,zSearch,1)
setappdata(handles.figure1,'lastClick',[newEstimatePoint(:,1:2)]);
hold(handles.axes1,'on')
plot(handles.axes1,cornersX,cornersY,'g')
hold(handles.axes1,'off');
drawnow




% newEstimatePoint=tpswarp3(affineFiducials,[],masterPoints,affineFiducials);
function reclick(hObject,eventdata)
handles=guidata(get(hObject,'Parent'));

newEstimatePoint=getappdata(handles.figure1,'lastClick');
windowSearch=str2double(get(handles.xySearch,'String'));
zSearch=str2double(get(handles.zSearch,'String'));
cornersX=newEstimatePoint(:,1)+[windowSearch windowSearch -windowSearch -windowSearch windowSearch];
cornersY=newEstimatePoint(:,2)+[windowSearch -windowSearch -windowSearch windowSearch windowSearch];
inputData(hObject,newEstimatePoint(:,1),newEstimatePoint(:,2),windowSearch,zSearch,1)
setappdata(handles.figure1,'lastClick',[newEstimatePoint(:,1:2)]);
hold(handles.axes1,'on')
plot(handles.axes1,cornersX,cornersY,'g')
hold(handles.axes1,'off');


function inputData(hObject,xselect,yselect,windowSearch,zSearch,manualFlag)
if nargin==5
    manualFlag=0;
end
manualFlag=double(manualFlag);

handles=guidata(get(hObject,'Parent'));
fiducialFile=get(handles.currentFiducialFile,'String');
iUser=get(handles.usersDropdown,'Value');
contents = cellstr(get(handles.usersDropdown,'String'));
iFrame=getappdata(handles.figure1,'currentFrame');
user=contents{iUser};
userFiducialsMaster=[];
userFiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
userFiducialsAll=userFiducialsAllUsers.(user);
plotIdx=getappdata(handles.figure1,'cursorTarget');

for iUser=1:length(contents)
    userFiducialsAllTemp=userFiducialsAllUsers.(contents{iUser});
    userFiducialsTemp=userFiducialsAllTemp{iFrame};
    if iUser==get(handles.usersDropdown,'Value')
        userFiducials=userFiducialsTemp;
        if size(userFiducialsTemp,1)>=plotIdx
            refPoint=userFiducialsTemp(plotIdx,:);
            userFiducialsTemp(plotIdx,:)=[];
            if isempty(refPoint{1})
                refPoint=[];
            end
        else
            refPoint=[];
        end
    end
    userFiducialsMaster=cat(1,userFiducialsMaster,userFiducialsTemp);
    
    
end

%     if isfield(fiducialData,'timeOffset');
% fiducialData = rmfield(fiducialData,'timeOffset');
%     end
%currentFiducials=(cellfun(@(x) cell2mat(x),fiducialData,'uniformOutput',0));

timeOffset=str2double(get(handles.timeOffset,'String'));
xRange=xlim(handles.axes1);
yRange=ylim(handles.axes1);
% if you click outside the image, the centroid will become nan
if xselect>xRange(1) && xselect< xRange(2) && yselect>yRange(1) && yselect<yRange(2);
    % turn on snapping later
    %     minD=pdist2([xselect,yselect],currentCentroids(:,1:2),'euclidean','smallest',1);
    %     pointIdx=find(minD==min(minD),1,'first');
    %     pointIdx=currentCentroids(pointIdx,3);
    %
    %baseImg=getappdata(0,'baseImg');
    filterSize=str2double(get(handles.filterSize,'String'));
    if ~isempty(filterSize)
        switch get(handles.filterOption,'Value')
            case 1
                if filterSize>0
                    gaussFilter=fspecial('gaussian', [10,10],filterSize);
                else
                    gaussFilter=1;
                end
                
            case 2
                gaussFilter=-fspecial('log', [10,10],filterSize);
                
        end
        
    else
        gaussFilter=fspecial('gaussian', [10,10],4);
        set(handles.filterSize,'String','4')
    end
    %baseImg=imfilter(baseImg,gaussFilter);
    %subSearch=baseImg(round(yselect)+(-windowSearch:windowSearch),round(xselect)+(-windowSearch:windowSearch));
    
    hiResIdx=getappdata(handles.figure1,'currentHiResIdx');
    frameIdx=getappdata(handles.figure1,'FrameIdx');
    
    %look in small volume around point
    if getappdata(handles.figure1,'tifFlag')
        temp=getappdata(handles.figure1,'imageStack');
        zRange=(-zSearch:zSearch)+hiResIdx;
        zRange(zRange<1)=1;
        maxVal=get(handles.zSlider,'Max');
        zRange(zRange>maxVal)=maxVal;
        worm=temp(:,:,zRange);
        
    else
    Fid=getappdata(handles.figure1,'fileID');

    
    if isempty(Fid)
        Fid=fopen([imFiles filesep 'sCMOS_Frames_U16_1024x1024.dat'] );
        setappdata(handles.figure1,'fileID',Fid);
    end
    R=getappdata(0,'registration');
    [row,col]=size(R.initialIm);
    status=fseek(Fid,2*(hiResIdx-zSearch)*row*col,-1);
    temp=fread(Fid,row*col*(2*zSearch+1),'uint16',0,'l');
    temp=(reshape(temp,row,col,(2*zSearch+1)));
    % fclose(Fid)
    %     temp=pixelIntensityCorrection(temp);
    %crop left and right regions
    rect1=R.rect1;
    rect2=R.rect2;
    t_concord=R.t_concord;
    Rsegment=R.Rsegment;
    padRegion=R.padRegion;
    
    
    if get(handles.channelSelect,'Value')==1
        worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    else
        worm=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
        worm=imwarp(worm,t_concord,'OutputView',Rsegment);
        
    end
    end
    worm=imfilter(worm,gaussFilter);
    subSearch=worm(round(yselect)+(-windowSearch:windowSearch),round(xselect)+(-windowSearch:windowSearch),:);
    %maxRegions=imregionalmax(subSearch);
    
    maxRegions=false(size(subSearch));
    for iSlice=1:size(subSearch,3);
        maxRegions(:,:,iSlice)=imregionalmax(subSearch(:,:,iSlice));
    end
    
    
    if ~any(maxRegions)
        maxRegions=(subSearch==max(subSearch));
    end
    
    [maxPosY,maxPosX,maxPosZ]=ind2sub(size(subSearch),find(maxRegions));
    
    
    zVoltages=getappdata(handles.figure1,'zVoltages');
    if mean(diff(frameIdx))>0
        subZVoltages=zVoltages(interp1(sort(frameIdx(:))+(1:length(frameIdx))'/100,...
            1:length(frameIdx),hiResIdx-timeOffset+(-zSearch:zSearch),'nearest',1));
    else
        subZVoltages=zVoltages(interp1(sort(frameIdx(:),'descend')+(length(frameIdx):-1:1)'/100,...
            1:length(frameIdx),hiResIdx-timeOffset+(-zSearch:zSearch),'nearest',1));
        
    end
    maxVals=subSearch(maxRegions);
    xselect=xselect-windowSearch+maxPosX-1;
    yselect=yselect-windowSearch+maxPosY-1;
    zselect=subZVoltages(maxPosZ);
    zSliceSelect=maxPosZ-1-zSearch+hiResIdx;
    if length(xselect)>1
        if size(userFiducialsMaster,2)<3
            currentXY=[0 0 0];
        else
            currentXY=cell2mat(userFiducialsMaster(:,1:4));
            
            if isempty(currentXY)
                currentXY=[0 0 0 0 0];
            end
            
        end
        
        %get only points that are further than 10 from all current points in the plane
%         dmat=pdist2([xselect,yselect], currentXY(:,1:2))<5;
%         eqmat=bsxfun(@minus ,zSliceSelect,currentXY(:,4)');
%         eqmat=abs(eqmat)<=1;
%         goodPoints=find(~any(dmat & eqmat,2));
%         if ~isempty(refPoint)
%             dmat2=pdist2([xselect(goodPoints),yselect(goodPoints)],...
%                 cell2mat(refPoint(1:2)))<1;
%             goodPoints=goodPoints(~dmat2);
%         end
%         if isempty(goodPoints)
%             goodPoints=1:length(maxVals);
%         end
%         
%         goodPoints=goodPoints((maxVals(goodPoints)==max(maxVals(goodPoints))));
        [~,goodPoints]=max(maxVals);
        xselect=xselect(goodPoints);
        yselect=yselect(goodPoints);
        zselect=zselect(goodPoints);
        maxPosZ=maxPosZ(goodPoints);
    end
    set(handles.zSlider,'Value',zselect)
    
    ctrlPnt=[xselect(1),yselect(1),zselect(1)];
    userFiducials{plotIdx,4}=getappdata(handles.figure1,'currentHiResIdx')-(zSearch+1)+maxPosZ(1);
    
    
    if isempty(userFiducials{plotIdx,1})
        clickPoints=pointUpdate(handles);
    else
        clickPoints=getappdata(handles.figure1,'points');
    end
    % you only get points if you are replacing an empty point
    
    userFiducials{plotIdx,1}=ctrlPnt(1);
    userFiducials{plotIdx,2}=ctrlPnt(2);
    userFiducials{plotIdx,3}=ctrlPnt(3);
    userFiducials{plotIdx,5}=manualFlag;
    
else
    clickPoints=getappdata(handles.figure1,'points');
    
    userFiducials{plotIdx,1}=[];
    userFiducials{plotIdx,2}=[];
    userFiducials{plotIdx,3}=[];
    userFiducials{plotIdx,4}=[];
    userFiducials{plotIdx,5}=[];
end

% if ~size(plotIdx
% end



userFiducialsAll{iFrame}=userFiducials;
set(handles.DisplayIdx,'data',userFiducials);
userFiducialsAllUsers.(user)=userFiducialsAll;

setappdata(handles.figure1,'fiducialsAllUsers',userFiducialsAllUsers);

setappdata(handles.figure1,'fiducials',userFiducialsAll);
fiducialPoints=userFiducialsAll;
setappdata(handles.figure1,'lastUser',user);

plotter(handles.slider1,'eventdata');

if getappdata(handles.figure1,'updateFlag')
    
    
    updataData2(handles)
end
setappdata(handles.figure1,'updateFlag',0);

if strcmp(get(handles.timer, 'Running'), 'on')
else               % If timer is stopped, reset its period only.
    start(handles.timer)
end



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

%any key press turns off cruise control
if get(handles.cruiseControl,'Value')
    % if strcmp(evnt.EventName,'KeyPress')
    set(handles.cruiseControl,'Value',0)
    % end
    return
end
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
elseif strcmp(evnt.Key,'space')
    cursorNeuronSelect(handles.slider1,evnt)
elseif strcmp(evnt.Key,'e');
    if get(handles.missingMode,'Value')
        nextMissing_Callback(hObject,evnt,handles)
    else
        goForward_Callback(handles.slider1,evnt,handles);
    end
    autoSelect(handles.slider1,evnt)
elseif strcmp(evnt.Key,'r')
    autoSelect(handles.slider1,evnt)
elseif strcmp(evnt.Key,'f');
    nextMissing_Callback(handles.slider1,evnt,handles)
elseif strcmp(evnt.Key,'z');
    previousMissing_Callback(handles.slider1,evnt,handles)
elseif strcmp(evnt.Key,'1')
    selectNeuron1_Callback(handles.selectNeuron1,evnt,handles);
    cursorNeuronSelect(handles.slider1,evnt)
elseif strcmp(evnt.Key,'2')
    selectNeuron2_Callback(handles.selectNeuron2,evnt,handles);
    cursorNeuronSelect(handles.slider1,evnt);
elseif strcmp(evnt.Key,'3')
    selectNeuron3_Callback(handles.selectNeuron3,evnt,handles);
    cursorNeuronSelect(handles.slider1,evnt)
elseif strcmp(evnt.Key,'4')
    selectNeuron4_Callback(handles.selectNeuron4,evnt,handles);
    cursorNeuronSelect(handles.slider1,evnt);
elseif strcmp(evnt.Key,'5')
    selectNeuron5_Callback(handles.selectNeuron5,evnt,handles);
    cursorNeuronSelect(handles.slider1,evnt);
elseif strcmp(evnt.Key,'6')
    selectNeuron6_Callback(handles.selectNeuron6,evnt,handles);
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
elseif strcmp(evnt.Key,'g')
    currentCursor=get(handles.selectNeuron1,'String');
    currentCursor=str2double(currentCursor)+1;
    set(handles.selectNeuron1,'String',num2str(currentCursor))
    set(handles.neuron1,'String',num2str(currentCursor))
    selectNeuron1_Callback(handles.selectNeuron1,evnt,handles);
    cursorNeuronSelect(handles.slider1,evnt)
end


%         elseif strcmp(evnt.Character,'h')
%             dispFeat=~dispFeat;
%             RefreshDisplayAndPlot;
%             disp('Hide/Show Features');
%         end
%
%Ignore the key stroke
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
imFiles=imFiles{1};
if ~isdir(imFiles)
    parent=fileparts(imFiles);
else
    parent=imFiles;
end
fiducialFile=uipickfiles('filterspec',parent);

fiducialFile=fiducialFile{1};
%fiducialData=load(fiducialFile);
timeOffset=load([fiducialFile filesep 'timeOffset']);
timeOffset=timeOffset.timeOffset;

set(handles.timeOffset,'String',timeOffset)

%setappdata(handles.figure1,'fiducialsData', fiducialData2)
set(handles.currentFiducialFile,'String',fiducialFile)
userList=dir([fiducialFile filesep '*.mat']); %get all user files
userList={userList.name};
userList=userList(~strcmp(userList,'timeOffset.mat'));
userList=cellfun(@(x) x(1:end-4),userList,'uniform',0);
handles.usersDropdown.String=userList;

pointUpdate(handles);
updataData2(handles);
if get(handles.savingFiducials,'Value')==0
    set(handles.savingFiducials,'Value',1);
    savingFiducials_Callback(handles.savingFiducials, eventdata, handles)
end

stop(handles.timer)

start(handles.timer);


function currentFiducialFile_Callback(hObject, eventdata, handles)
% hObject    handle to currentFiducialFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentFiducialFile as text
%        str2double(get(hObject,'String')) returns contents of currentFiducialFile as a double


% --- Executes during object creation, after setting all properties.
function currentFiducialFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentFiducialFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previousAnnotated.
function previousAnnotated_Callback(hObject, eventdata, handles)
% hObject    handle to previousAnnotated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iUser=get(handles.usersDropdown,'Value');
contents = cellstr(get(handles.usersDropdown,'String'));
user=contents{iUser};
userFiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
fiducialPoints=userFiducialsAllUsers.(user);
currentFrame=getappdata(handles.figure1,'currentFrame');
plotIdx=getappdata(handles.figure1,'cursorTarget');
ciLimit= str2double(get(handles.confidenceLimit,'String'));

try
    W=cellfun(@(x) (x{plotIdx,5}),fiducialPoints,'uniformOutput',0);
    missing=cellfun(@(x) isempty(x{plotIdx,1}),fiducialPoints);
    W(missing)={0};
    W=cell2mat(W);
    missing=find(W>=ciLimit);
catch
    missing=(cellfun(@(x) size(x,1)<plotIdx,fiducialPoints));
    Wpresent=cellfun(@(x) (x{plotIdx,5}),fiducialPoints(~missing),'uniform',0);
    Wpresent2=cellfun(@(x) ~isempty(x),Wpresent);
    Wpresent3=double(Wpresent2);
Wpresent3(Wpresent2)=cell2mat(Wpresent);
    W=zeros(size(missing));
    W(~missing)=Wpresent3;
    missing=find(W>=ciLimit);
end
nextFrame=missing((missing<currentFrame));
if isempty(nextFrame)
    return
end
nextFrame=nextFrame(end);
moveFrame(hObject,eventdata,handles,nextFrame-currentFrame)



% --- Executes on button press in nextAnnotated.
function nextAnnotated_Callback(hObject, eventdata, handles)
% hObject    handle to nextAnnotated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


iUser=get(handles.usersDropdown,'Value');
contents = cellstr(get(handles.usersDropdown,'String'));
user=contents{iUser};
userFiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
fiducialPoints=userFiducialsAllUsers.(user);
currentFrame=getappdata(handles.figure1,'currentFrame');
plotIdx=getappdata(handles.figure1,'cursorTarget');
ciLimit= str2double(get(handles.confidenceLimit,'String'));

try
    W=cellfun(@(x) (x{plotIdx,5}),fiducialPoints,'uniformOutput',0);
    missing=cellfun(@(x) isempty(x{plotIdx,1}),fiducialPoints);
    W(missing)={0};
    W=cell2mat(W);
    missing=find(W>=ciLimit);
catch
    missing=(cellfun(@(x) size(x,1)<plotIdx,fiducialPoints));
    Wpresent=cellfun(@(x) (x{plotIdx,5}),fiducialPoints(~missing),'uniform',0);
    Wpresent2=cellfun(@(x) ~isempty(x),Wpresent);
    Wpresent3=double(Wpresent2);
Wpresent3(Wpresent2)=cell2mat(Wpresent);
    W=zeros(size(missing));
    W(~missing)=Wpresent3;
    missing=find(W>=ciLimit);
end
nextFrame=missing((missing>currentFrame));
if isempty(nextFrame)
    return
end
nextFrame=nextFrame(1);
moveFrame(hObject,eventdata,handles,nextFrame-currentFrame)


% --- Executes on button press in previousMissing.
function previousMissing_Callback(hObject, eventdata, handles)
% hObject    handle to previousMissing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iUser=get(handles.usersDropdown,'Value');
contents = cellstr(get(handles.usersDropdown,'String'));
user=contents{iUser};
userFiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
fiducialPoints=userFiducialsAllUsers.(user);
currentFrame=getappdata(handles.figure1,'currentFrame');
plotIdx=getappdata(handles.figure1,'cursorTarget');
ciLimit= str2double(get(handles.confidenceLimit,'String'));

try
    W=cellfun(@(x) (x{plotIdx,5}),fiducialPoints,'uniformOutput',0);
    missing=cellfun(@(x) isempty(x{plotIdx,1}),fiducialPoints);
    W(missing)={0};
    W=cell2mat(W);
    missing=find(W<ciLimit);
catch
    missing=(cellfun(@(x) size(x,1)<plotIdx,fiducialPoints));
    Wpresent=cellfun(@(x) (x{plotIdx,5}),fiducialPoints(~missing),'uniform',0);
    Wpresent2=cellfun(@(x) ~isempty(x),Wpresent);
Wpresent2(Wpresent2)=cell2mat(Wpresent);
    W=zeros(size(missing));
    W(~missing)=Wpresent2;
    missing=find(W<ciLimit);
end
nextFrame=missing((missing<currentFrame));
if isempty(nextFrame)
    return
end
nextFrame=nextFrame(end);
moveFrame(hObject,eventdata,handles,nextFrame-currentFrame)


% --- Executes on button press in nextMissing.
function nextMissing_Callback(hObject, eventdata, handles)
% hObject    handle to nextMissing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iUser=get(handles.usersDropdown,'Value');
contents = cellstr(get(handles.usersDropdown,'String'));
user=contents{iUser};
userFiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
fiducialPoints=userFiducialsAllUsers.(user);
currentFrame=getappdata(handles.figure1,'currentFrame');
plotIdx=getappdata(handles.figure1,'cursorTarget');
ciLimit= str2double(get(handles.confidenceLimit,'String'));

try
    W=cellfun(@(x) (x{plotIdx,5}),fiducialPoints,'uniformOutput',0);
    missing=cellfun(@(x) isempty(x{plotIdx,1}),fiducialPoints);
    W(missing)={0};
    W=cell2mat(W);
    missing=find(W<ciLimit);
catch
    missing=(cellfun(@(x) size(x,1)<plotIdx,fiducialPoints));
    Wpresent=cellfun(@(x) (x{plotIdx,5}),fiducialPoints(~missing),'uniform',0);
    Wpresent2=(cellfun(@(x) ~isempty(x),Wpresent));
    Wpresent3=double(Wpresent2);
Wpresent3(Wpresent2)=cell2mat(Wpresent);
    W=zeros(size(missing));
    W(~missing)=Wpresent3;
    missing=find(W<ciLimit);
end
nextFrame=missing((missing>currentFrame));
if isempty(nextFrame)
    return
end
nextFrame=nextFrame(1);
moveFrame(hObject,eventdata,handles,nextFrame-currentFrame)

% set(handles.slider1,'Value',min(nextFrame,get(handles.slider1,'max')));
% plotter(handles.slider1,eventdata);



% --- Executes on button press in savingFiducials.
function savingFiducials_Callback(hObject, eventdata, handles)
% hObject    handle to savingFiducials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'value')
    set(hObject,'String','SAVING','fontsize', 26, 'BackgroundColor',[1,0,0] )
else
    set(hObject,'String','save fiducials','fontsize', 12, 'BackgroundColor',[1,1,1] )
end


% --- Executes on button press in Continuous.
function Continuous_Callback(hObject, eventdata, handles)
% hObject    handle to Continuous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Continuous


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.figure1,'SelectionType'),'alt')
    cursorNeuronSelect(handles.slider1,eventdata);
    
end


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



%updates the number of points, possibly playing fun sounds
function points=pointUpdate(handles)
currentClick=datevec(now);
history=getappdata(handles.figure1,'history');
points=getappdata(handles.figure1,'points');
multiplier=1;


if size(history,1)>30
    history=cat(1,history(2:end,:),currentClick);
else
    history=cat(1,history,currentClick);
end
setappdata(handles.figure1,'history',history)
if size(history,1)>10
    
    timeIntervalHistory=diff(history,[],1);
    timeIntervalHistory=timeIntervalHistory*[0 0 3600*24 3600 60 1]';
    
    shortTimeHistory=mean(timeIntervalHistory(end-9:end));
    if size(history,1)>20
        longTimeHistory=mean(timeIntervalHistory);
    else
        longTimeHistory=inf;
    end
    %multiplier can be up to 10.
    if shortTimeHistory<10
        multiplier=multiplier+1;
    end
    if shortTimeHistory<5
        multiplier=multiplier+3;
    end
    
    if longTimeHistory<10
        multiplier=multiplier+4;
    end
    if longTimeHistory<5
        multiplier=multiplier+6;
    end
    
    if longTimeHistory<3.5
        multiplier=multiplier+5;
    end
    %[shortTimeHistory multiplier]
    
end
fireCombo=getappdata(handles.figure1,'fireCombo');
switch multiplier
    case {2 5}
        Etext='Keep it up!';
        fireCombo=0;
    case {6 9}
        Etext='You''re on a roll!';
        fireCombo=0;
    case {15}
        Etext='Amazing!!';
        fireCombo=0;
    case {20}
        Etext='You''re on FIRE!!!';
        if fireCombo==0
            soundpath=which('Worm3DFiducialPickerNetwork');
            soundpath=fileparts(soundpath);
            soundpath=[soundpath filesep 'onFire.mp3'];
            [data,fs]=audioread(soundpath);
            sound(data,fs);
        end
        
        fireCombo=fireCombo+1/20;
        if floor(fireCombo)>0
            Etext=[Etext '   Fire Combo x' num2str(floor(fireCombo))];
            
        end
        
        if abs(fireCombo-round(fireCombo))<.01
            soundpath=which('Worm3DFiducialPickerNetwork');
            soundpath=fileparts(soundpath);
            
            switch mod(round(fireCombo),34)
                
                case 3
                    soundpath=[soundpath filesep 'tripleCombo.mp3'];
                case 6
                    soundpath=[soundpath filesep 'superCombo.mp3'];
                case 9
                    soundpath=[soundpath filesep 'hyperCombo.mp3'];
                case 12
                    soundpath=[soundpath filesep 'brutalCombo.mp3'];
                case 15
                    soundpath=[soundpath filesep 'masterCombo.mp3'];
                case 18
                    soundpath=[soundpath filesep 'awesomeCombo.mp3'];
                case 21
                    soundpath=[soundpath filesep 'blasterCombo.mp3'];
                case 24
                    soundpath=[soundpath filesep 'monsterCombo.mp3'];
                case 27
                    soundpath=[soundpath filesep 'kingCombo.mp3'];
                case 30
                    soundpath=[soundpath filesep 'killerCombo.mp3'];
                case 33
                    soundpath=[soundpath filesep 'ultraCombo.mp3'];
                otherwise
                    if rand>.7
                        
                        soundpath=[soundpath filesep 'Boomshakalaka.wav'];
                    elseif rand>.5
                        soundpath=[soundpath filesep 'kaboom.mp3'];
                    elseif rand>.3
                        soundpath=[soundpath filesep 'boomshakalaka2.mp3'];
                    else
                        soundpath=[soundpath filesep 'onFire.mp3'];
                    end
                    
            end
            
            [data,fs]=audioread(soundpath);
            sound(data,fs);
            display('Boomshakalaka!')
            
        end
        
    otherwise
        Etext='Click the Neurons!';
        fireCombo=0;
end

if fireCombo<getappdata(handles.figure1,'fireCombo');
    Etext=['COMBO BREAKER'];
    soundpath=which('Worm3DFiducialPickerNetwork');
    soundpath=fileparts(soundpath);
    [data,fs]=audioread([soundpath filesep 'combobreaker.mp3']);
    sound(data,fs);
    display('COMBO BREAKER');
end
setappdata(handles.figure1,'fireCombo',fireCombo(1));

totalNeuronsClicked=getappdata(handles.figure1,'totalNeuronsClicked');

Etext=[ num2str(totalNeuronsClicked) ' Total Neurons.    ' Etext];

set(handles.encouragingText,'string',Etext);

points=points+multiplier;
setappdata(handles.figure1,'points',points);
level=floor(log(points/10000)/log(1.3))+1;
level=max(level,0);
if isempty(getappdata(handles.figure1,'level'));
    setappdata(handles.figure1,'level',level);
    oldLevel=level;
else
    oldLevel=getappdata(handles.figure1,'level');
end

if level>oldLevel;
    try
        pointsAll=getappdata(handles.figure1,'pointsAll');
        iUser=get(handles.usersDropdown,'Value');
        [~,ia]=sort(pointsAll,'descend');
        
        display(['Congratulations! You have just reached level: ' ...
            num2str(level) '. Your rank is ' num2str(find(ia==iUser))]);
    catch
    end
    
end
setappdata(handles.figure1,'level',level);


points2nextLevel=round(10000*(1.3)^(level)-points);

set(handles.Points,'String',['Level:' num2str(level) ...
    '| Points to next level: ' num2str(points2nextLevel)]);



function xySearch_Callback(hObject, eventdata, handles)
% hObject    handle to xySearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xySearch as text
%        str2double(get(hObject,'String')) returns contents of xySearch as a double


% --- Executes during object creation, after setting all properties.
function xySearch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xySearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zSearch_Callback(hObject, eventdata, handles)
% hObject    handle to zSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zSearch as text
%        str2double(get(hObject,'String')) returns contents of zSearch as a double


% --- Executes during object creation, after setting all properties.
function zSearch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function refIdx_Callback(hObject, eventdata, handles)
% hObject    handle to refIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refIdx as text
%        str2double(get(hObject,'String')) returns contents of refIdx as a double


% --- Executes during object creation, after setting all properties.
function refIdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filterSize_Callback(hObject, eventdata, handles)
% hObject    handle to filterSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterSize as text
%        str2double(get(hObject,'String')) returns contents of filterSize as a double


% --- Executes during object creation, after setting all properties.
function filterSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in filterOption.
function filterOption_Callback(hObject, eventdata, handles)
% hObject    handle to filterOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filterOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filterOption


% --- Executes during object creation, after setting all properties.
function filterOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function neuron1_Callback(hObject, eventdata, handles)
% hObject    handle to neuron1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron1 as text
%        str2double(get(hObject,'String')) returns contents of neuron1 as a double
currentNeuron=str2double(get(hObject,'String'));
currentNeuronString=num2str(currentNeuron);
set(handles.selectNeuron1,'String',currentNeuronString);
setappdata(handles.figure1,'cursorTarget',currentNeuron);
buttonClear(handles)
set(handles.selectNeuron1,'BackgroundColor',[1 0 0 ]);



function neuron2_Callback(hObject, eventdata, handles)
% hObject    handle to neuron2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron2 as text
%        str2double(get(hObject,'String')) returns contents of neuron2 as a double
currentNeuron=str2double(get(hObject,'String'));
currentNeuronString=num2str(currentNeuron);
set(handles.selectNeuron2,'String',currentNeuronString);
setappdata(handles.figure1,'cursorTarget',currentNeuron);
buttonClear(handles)
set(handles.selectNeuron2,'BackgroundColor',[1 0 0 ]);



function neuron3_Callback(hObject, eventdata, handles)
% hObject    handle to neuron3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron3 as text
%        str2double(get(hObject,'String')) returns contents of neuron3 as a double
currentNeuron=str2double(get(hObject,'String'));
currentNeuronString=num2str(currentNeuron);
set(handles.selectNeuron3,'String',currentNeuronString);
setappdata(handles.figure1,'cursorTarget',currentNeuron);
buttonClear(handles)
set(handles.selectNeuron3,'BackgroundColor',[1 0 0 ]);



function neuron4_Callback(hObject, eventdata, handles)
% hObject    handle to neuron4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron4 as text
%        str2double(get(hObject,'String')) returns contents of neuron4 as a double
currentNeuron=str2double(get(hObject,'String'));
currentNeuronString=num2str(currentNeuron);
set(handles.selectNeuron4,'String',currentNeuronString);
setappdata(handles.figure1,'cursorTarget',currentNeuron);
buttonClear(handles)
set(handles.selectNeuron4,'BackgroundColor',[1 0 0 ]);



function neuron5_Callback(hObject, eventdata, handles)
% hObject    handle to neuron5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron5 as text
%        str2double(get(hObject,'String')) returns contents of neuron5 as a double
currentNeuron=str2double(get(hObject,'String'));
currentNeuronString=num2str(currentNeuron);
set(handles.selectNeuron5,'String',currentNeuronString);
setappdata(handles.figure1,'cursorTarget',currentNeuron);
buttonClear(handles)
set(handles.selectNeuron5,'BackgroundColor',[1 0 0 ]);




function neuron6_Callback(hObject, eventdata, handles)
% hObject    handle to neuron6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron6 as text
%        str2double(get(hObject,'String')) returns contents of neuron6 as a double
currentNeuron=str2double(get(hObject,'String'));
currentNeuronString=num2str(currentNeuron);
set(handles.selectNeuron6,'String',currentNeuronString);
setappdata(handles.figure1,'cursorTarget',currentNeuron);
buttonClear(handles)
set(handles.selectNeuron6,'BackgroundColor',[1 0 0 ]);


function buttonClear(handles)
%sets all buttons for neuron selection to 0 state
set(handles.selectNeuron1,'Value',0);
set(handles.selectNeuron2,'Value',0);
set(handles.selectNeuron3,'Value',0);
set(handles.selectNeuron4,'Value',0);
set(handles.selectNeuron5,'Value',0);
set(handles.selectNeuron6,'Value',0);
set(handles.selectNeuron1,'BackgroundColor',[.9 .9 .9 ]);
set(handles.selectNeuron2,'BackgroundColor',[.9 .9 .9 ]);
set(handles.selectNeuron3,'BackgroundColor',[.9 .9 .9 ]);
set(handles.selectNeuron4,'BackgroundColor',[.9 .9 .9 ]);
set(handles.selectNeuron5,'BackgroundColor',[.9 .9 .9 ]);
set(handles.selectNeuron6,'BackgroundColor',[.9 .9 .9 ]);




% --- Executes during object creation, after setting all properties.
function neuron5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function neuron6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function neuron1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function neuron2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function neuron3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function neuron4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in usersDropdown.
function usersDropdown_Callback(hObject, eventdata, handles)
% hObject    handle to usersDropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns usersDropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from usersDropdown
contents = cellstr(get(hObject,'String'));
user=contents{get(hObject,'Value')};
fiducialFile=get(handles.currentFiducialFile,'String');
updataData2(handles);


% --- Executes during object creation, after setting all properties.
function usersDropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usersDropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pop=[fileparts(which('Worm3dFiducialPickerNetwork.m')) filesep 'worm3dFiducialPickerNetworkHelp.txt'];
popupmessage(pop,'Worm3dFiducialPickerNetwork Controls')

function updateData(hObject,eventdata,hfigure)
%periodically update data from selected fiducial folder
handles = guidata(hfigure);

setappdata(handles.figure1,'updateFlag',1);




function updataData2(handles)
contents = cellstr(get(handles.usersDropdown,'String'));
fiducialFileName=get(handles.currentFiducialFile,'String');
successFlag=false;
tryCounter=0;
ME=[];
display('Refreshing Data');
currentUser=contents(get(handles.usersDropdown,'Value'));
userList=dir([fiducialFileName filesep '*.mat']); %get all user files
userList={userList.name};
userList=userList(~strcmp(userList,'timeOffset.mat'));
userList=cellfun(@(x) x(1:end-4),userList,'uniform',0);
handles.usersDropdown.String=userList;

pointUpdate(handles);
%updataData2(handles);

%save current user
if get(handles.savingFiducials,'Value') && ~isempty(getappdata(handles.figure1,'lastUser'))
    clickPoints=getappdata(handles.figure1,'points');
    userFiducialsAll=getappdata(handles.figure1,'fiducialsAllUsers');
    user=getappdata(handles.figure1,'lastUser');
    fiducialPoints=userFiducialsAll.(user);
    fiducialFile=get(handles.currentFiducialFile,'String');
    timeOffset=str2double(get(handles.timeOffset,'String'));
    save([fiducialFile filesep user],'fiducialPoints','clickPoints','-v6')
    save([fiducialFile filesep 'timeOffset'],'timeOffset','-v6')
end



while ~successFlag
    tryCounter=tryCounter+1;

    try
        % loop through and load users from file
        for iUser=1:length(userList);
            
            
            user=userList{iUser};
            if exist([fiducialFileName filesep user '.mat'],'file')
                fiducialsAll=load([fiducialFileName filesep user]);
                if isfield(fiducialsAll,'clickPoints')
                points=fiducialsAll.clickPoints;
                else
                    points=1;
                end
                fiducialsAll=fiducialsAll.fiducialPoints;
                 fiducialsAll=fiducialsAll(:);
%             else
%                 hiResData=getappdata(handles.figure1,'hiResData');
%                 fiducialPoints=cell(6,5);
%                 fiducialPoints=repmat({fiducialPoints},max(hiResData.stackIdx),1);
%                 setappdata(handles.figure1,'fiducials',fiducialPoints);
%                 clickPoints=0;
%                 points=clickPoints;
%                 save([fiducialFileName filesep user],'fiducialPoints','clickPoints','-v6');
%                 fiducialsAll=fiducialPoints;
%                 
            end
            
            if strcmp(currentUser,user)
                setappdata(handles.figure1,'points',points);
                set(handles.usersDropdown,'Value',iUser)
            end
            pointsAll(iUser)=points;
            
            fiducialsAllUsers.(user)=fiducialsAll;
            successFlag=true;
        end
        setappdata(handles.figure1,'pointsAll',pointsAll);
        setappdata(handles.figure1,'fiducialsAllUsers',fiducialsAllUsers);
        nFidClicked=@(Y) max(cell2mat(cellfun(@(x) size(cell2mat(x(:,1)),1),Y,'uniformOutput',0)));
            try
        totalNeuronsClicked=structfun(nFidClicked, fiducialsAllUsers,'uniform',0);
        totalNeuronsClicked=sum(cell2mat(struct2cell(totalNeuronsClicked)));
            catch
                totalNeuronsClicked=0;
            end
        setappdata(handles.figure1,'totalNeuronsClicked',totalNeuronsClicked);
        
    catch ME
        successFlag=false;
        if tryCounter>5
            display('ERROR, Fiducial Loading fail')
            rethrow(ME)
        else
            display('Loading Conflict, trying again');
        end
        
        
    end
end
display('Refresh Complete');



function refreshRate_Callback(hObject, eventdata, handles)
% hObject    handle to refreshRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refreshRate as text
%        str2double(get(hObject,'String')) returns contents of refreshRate as a double
period= str2double(get(hObject,'String')) ;
if strcmp(get(handles.timer, 'Running'), 'on')
    stop(handles.timer);
    set(handles.timer,'Period',period)
    start(handles.timer)
else               % If timer is stopped, reset its period only.
    set(handles.timer,'Period',period)
    start(handles.timer)
end

% --- Executes during object creation, after setting all properties.
function refreshRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refreshRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updateButton.
function updateButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updataData2(handles)


% --- Executes on button press in cruiseControl.
function cruiseControl_Callback(hObject, eventdata, handles)
% hObject    handle to cruiseControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentFrame=getappdata(handles.figure1,'currentFrame');
setappdata(handles.figure1,'cruiseStartFrame',currentFrame);

%while get(hObject,'Value')
while get(hObject,'Value')
    if get(handles.missingMode,'Value')
        nextMissing_Callback(hObject,eventdata,handles)
    else
        goForward_Callback(handles.slider1,eventdata,handles);
    end
    autoSelect(handles.slider1,eventdata)
end

setappdata(handles.figure1,'cruiseStartFrame',[]);


% --- Executes on button press in errorCheck.
function errorCheck_Callback(hObject, eventdata, handles)
% hObject    handle to errorCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fiducialFileName=get(handles.currentFiducialFile,'String');

pop=[fiducialFileName filesep 'ErrorNotes.txt'];
if exist(pop)==2;
    
    popupmessage(pop,'Error Notes')
else
    fid=fopen(pop,'w');
    fclose(fid);
end
open(pop);


% --- Executes on button press in missingMode.
function missingMode_Callback(hObject, eventdata, handles)
% hObject    handle to missingMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of missingMode



function confidenceLimit_Callback(hObject, eventdata, handles)
% hObject    handle to confidenceLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of confidenceLimit as text
%        str2double(get(hObject,'String')) returns contents of confidenceLimit as a double


% --- Executes during object creation, after setting all properties.
function confidenceLimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to confidenceLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in renameDataset.
function renameDataset_Callback(hObject, eventdata, handles)
% hObject    handle to renameDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userList={'Fred' ,'Ashley','David','Mochi',...
    'Sagar','George','Jeff','Andy','Kevin','Lukas', 'Jose','Josh'};
contents = cellstr(get(handles.usersDropdown,'String'));
fiducialFileName=get(handles.currentFiducialFile,'String');
successFlag=false;
tryCounter=0;
ME=[];
display('Refreshing Data');
    user=get(handles.usersDropdown,'Value');
    stringList=handles.usersDropdown.String;


% if get(handles.savingFiducials,'Value') && ~isempty(getappdata(handles.figure1,'lastUser'))
%     clickPoints=getappdata(handles.figure1,'points');
%     userFiducialsAll=getappdata(handles.figure1,'fiducialsAllUsers');
%     fiducialPoints=userFiducialsAll.(user);
%     fiducialFile=get(handles.currentFiducialFile,'String');
%     timeOffset=str2double(get(handles.timeOffset,'String'));
%     save([fiducialFile filesep user],'fiducialPoints','clickPoints','-v6')
%     save([fiducialFile filesep 'timeOffset'],'timeOffset','-v6')
% end
choice = menu('Choose a user',userList);
destination=[fiducialFileName filesep userList{choice}];
nNumber=length(dir([destination '*']));
newUser=[userList{choice}  num2str(nNumber+1)]; 
destination=[destination num2str(nNumber+1) '.mat'];
source=[fiducialFileName filesep contents{user} '.mat'];
stringList(strcmp(stringList,contents{user}))={newUser};
movefile(source,destination)
handles.usersDropdown.String=stringList;
updataData2(handles)


% --- Executes on button press in compareOnlyPrevious.
function compareOnlyPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to compareOnlyPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compareOnlyPrevious


% --- Executes on button press in deleteNeuron.
function deleteNeuron_Callback(hObject, eventdata, handles)
% hObject    handle to deleteNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
killIdx=inputdlg('Which Neuron should be deleted?');
killIdx=str2double(killIdx);
iUser=get(handles.usersDropdown,'Value');
contents = cellstr(get(handles.usersDropdown,'String'));
user=contents{iUser};
userFiducialsAllUsers=getappdata(handles.figure1,'fiducialsAllUsers');
fiducialPoints=userFiducialsAllUsers.(user);

fiducialPoints=deleteNeuron(fiducialPoints,killIdx);
userFiducialsAllUsers.(user)=fiducialPoints;
setappdata(handles.figure1,'fiducialsAllUsers',userFiducialsAllUsers);
setappdata(handles.figure1,'lastUser',user)
updataData2(handles);
