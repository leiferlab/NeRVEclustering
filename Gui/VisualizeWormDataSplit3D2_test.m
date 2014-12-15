function varargout = VisualizeWormDataSplit3D2_test(varargin)
% VISUALIZEWORMDATASPLIT3D2_TEST MATLAB code for VisualizeWormDataSplit3D2_test.fig
%      VISUALIZEWORMDATASPLIT3D2_TEST, by itself, creates a new VISUALIZEWORMDATASPLIT3D2_TEST or raises the existing
%      singleton*.
%
%      H = VISUALIZEWORMDATASPLIT3D2_TEST returns the handle to a new VISUALIZEWORMDATASPLIT3D2_TEST or the handle to
%      the existing singleton*.
%
%      VISUALIZEWORMDATASPLIT3D2_TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZEWORMDATASPLIT3D2_TEST.M with the given input arguments.
%
%      VISUALIZEWORMDATASPLIT3D2_TEST('Property','Value',...) creates a new VISUALIZEWORMDATASPLIT3D2_TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisualizeWormDataSplit3D2_test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisualizeWormDataSplit3D2_test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisualizeWormDataSplit3D2_test

% Last Modified by GUIDE v2.5 29-Oct-2014 00:59:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VisualizeWormDataSplit3D2_test_OpeningFcn, ...
    'gui_OutputFcn',  @VisualizeWormDataSplit3D2_test_OutputFcn, ...
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


% --- Executes just before VisualizeWormDataSplit3D2_test is made visible.
function VisualizeWormDataSplit3D2_test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualizeWormDataSplit3D2_test (see VARARGIN)

% Choose default command line output for VisualizeWormDataSplit3D2_test

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

% UIWAIT makes VisualizeWormDataSplit3D2_test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualizeWormDataSplit3D2_test_OutputFcn(hObject, eventdata, handles)
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
imFolder=getappdata(0,'imFolder');

%load mat file folder
if isempty(imFolder)
    imFolder = uigetdir([],'Select MatFile Folder');
    setappdata(0,'imFolder',imFolder);
else
    try
        imFolder = uigetdir(imFolder,'Select MatFile Folder');
        setappdata(0,'imFolder',imFolder);
    catch
        imFolder = uigetdir([],'Select MatFile Folder');
        setappdata(0,'imFolder',imFolder);
    end
end

%select image folder
display('Select image folder, select one folder for split image or select red then green folder');
rawImFolder=uipickfiles('FilterSpec', fileparts(imFolder));

if length(rawImFolder)==1
%load registration file if needed for split
[rpath,parent]=uigetfile('Y:\CommunalCode\3dbrain\','Select Registration File');
registration=load([parent filesep rpath]);
setappdata(0,'registration',registration);
end

matFiles=dir([imFolder filesep '*.mat']);
setappdata(0,'matFiles',matFiles);


for iFolder=1:length(rawImFolder)
imFiles{iFolder}=dir([rawImFolder{iFolder} filesep '*.tif']);
end

setappdata(0,'imFiles',imFiles);
setappdata(0,'rawImFolder',rawImFolder);

%setting slider parameters
set(handles.slider1,'Min',1)
if isempty(matFiles)
    set(handles.slider1,'Max',length(imFiles{1}));
else
    set(handles.slider1,'Max',length(imFiles{1}));
end

set(handles.slider1,'Value',1)

if ~exist([imFolder filesep 'trackOutput.mat'],'file')
    runTrack_Callback(hObject, eventdata, handles)
    
end
load([imFolder filesep 'trackOutput']);



setappdata(handles.figure1,'trackOutput',trackOutput);
setappdata(handles.figure1,'currentFrame',1);
set(handles.slider1,'value',1);
set(handles.slider2,'min',1);
set(handles.slider2,'max',max(trackOutput(:,end)));
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

%plots current frame in image window and plot window


handles=guidata(get(hObject,'Parent'));
timeStep=str2double(get(handles.timeStep,'string'));
smoothWindow=str2double(get(handles.smoothingWindow,'String'));
startTime=str2double(get(handles.startTime,'String'));
normalizeFlag=get(handles.normalizeButton,'value');
imFolder=getappdata(0,'imFolder');
iImage=round(get(handles.slider1,'Value'));


%load track data, and lists of mat and tif files
trackData=getappdata(handles.figure1,'trackOutput');

matFiles=getappdata(0,'matFiles');
imFiles=getappdata(0,'imFiles');
rawImFolder=getappdata(0,'rawImFolder');
try
    matName= matFiles(iImage).name;
    wormData=load([imFolder filesep matName]);
catch
    wormData=[];
    matName=['stack' num2str(iImage,'%04d') 'data'];
end

flag3d=0;



%proceed only if mat file has wormMask
if isfield(wormData,'wormMask')
    wormMask=wormData.wormMask;
else
    wormMask=ones(2,2,2); %if not fitted, assume matrix for now
    
end

if length(imFiles)==1
    R=getappdata(0,'registration');
    
    imFiles=imFiles{1};
    %load image and correct pixels
    temp=double(imread([rawImFolder filesep imFiles(iImage).name],'tif'));
    %     temp=pixelIntensityCorrection(temp);
    %crop left and right regions
    rect1=R.rect1;
    rect2=R.rect2;
    t_concord=R.t_concord;
    Rsegment=R.Rsegment;
    padRegion=R.padRegion;
    temp_activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
    worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
    
    %align 2 halves
    temp_activity=imwarp(temp_activity,t_concord,'OutputView',Rsegment);
    temp_activity(padRegion)=median(temp_activity(~padRegion));
    activity=(temp_activity);%bpass_jn(temp_activity,1,[40,40]);
    hold(handles.axes1,'off')
    %clear current axes
    arrayfun(@(x) delete(x),get(handles.axes1,'children'))
    
    %chose which image to display, red or green
    switch get(handles.channelSelect,'value')
        case 1
            baseImg=worm; %red
        case 2
            baseImg=activity; %green
    end
    
else %if imfolder has 2, then just lead the correct iage folder
    if ismatrix(wormMask)
        imFilesIdx=find(cellfun(@(x) ~isempty(x),strfind({imFiles{1}.name},matFiles(iImage).name(6:10))'));
        flag3d=0;
        
        switch get(handles.channelSelect,'value')
            case 1
                baseImg= double(imread([rawImFolder{1} filesep imFiles{1}(imFilesIdx).name]...
                    ,'tif')); %red
            case 2
                baseImg= double(imread([rawImFolder{2} filesep imFiles{2}(imFilesIdx).name]...
                    ,'tif')); %green
                
            otherwise
                baseImg=0;
        end
        imName=imFiles{2}(imFilesIdx).name;
        zSlice=1;
    else
        imFilesIdxStr=matName(strfind(matName,'k')+1:strfind(matName,'d')-1);
        imFilesIdx=str2num(imFilesIdxStr);
        try
        stackName=(['image' imFilesIdxStr]);
        metaData=load([rawImFolder{3} filesep stackName]);
        catch
        stackName=(['image0' imFilesIdxStr]);
        metaData=load([rawImFolder{3} filesep stackName]);       
        end
        
        metaData=metaData.metaData;
        if isfield(metaData,'metaData');
            metaData=metaData.metaData;
        end
        
        zPos=metaData.zVoltage;
        if isfield(metaData,'midPlane')
            midPlane=metaData.midPlane;
        else
            midPlane=round(length(zPos)/2);
        end
        
        
        zPos=zPos-interp1(zPos,midPlane);
        set(handles.zSlider,'Max',max(zPos))
        set(handles.zSlider,'Min',min(zPos))
        zPosV=(get(handles.zSlider,'Value'));
        if zPosV>max(zPos);
            zPosV=max(zPos);
            set(handles.zSlider,'Value',zPosV);
        end
        if zPosV<min(zPos)
            zPosV=min(zPos);
            set(handles.zSlider,'Value',zPosV);
        end
      %  zPos=sort(zPos);
        zPos=zPos+(1:length(zPos))'*.00001;
        zSlice=interp1((zPos),1:length(zPos),zPosV,'nearest');
        zSlice=max(zSlice,1);
        flag3d=1;
        imName=['image' num2str(metaData.iFrame(zSlice),'%3.5d')];
        try
            
                        labelMask=wormMask(:,:,zSlice);
        catch
        end
        
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
                centroids=wormData.centroids;
                centroids(:,3)=wormData.zPlaneIdx;
        end
        
    end
    
end
setappdata(0,'baseImg',baseImg)
%     figure
%     imagesc(smooth2a(baseImg,20,20)>5);
%   set(handles.FrameIdx,'string',[num2str(iImage*timeStep,'%6.2f') 's' ...
%         '  ' num2str(zPosV)]);
     set(handles.FrameIdx,'string',stackName);

if ~isempty(trackData);
if ~isempty(baseImg(:))
    if get(handles.channelSelect,'value')~=3
        baseImg=pedistalSubtract(baseImg);
        baseImg(isnan(baseImg))=0;
        
        setappdata(handles.figure1,'baseImg',baseImg);
        %scale dynamic range
        newContrast=getappdata(handles.figure1,'newContrast');
        if isempty(newContrast)
            newContrast=[min(baseImg(:)),max(baseImg(:))];
        end
        baseImg(baseImg<newContrast(1)) = newContrast(1);
        baseImg(baseImg>newContrast(2)) = newContrast(2);
        baseImg = (baseImg-newContrast(1))./diff(newContrast);
    end
        hold(handles.axes1,'off')

    ax1=imagesc(baseImg,'parent',handles.axes1);
    
    
    plot(handles.axes3,smooth(max(smooth2a(baseImg,10,10),[],2),30));
    hold(handles.axes1,'on')
    
    %scatter centroids and label
    
    tracks=trackData(trackData(:,end-1)==iImage,:);
    
    %only plot centroids that appear in the plane of interest
    tracks=tracks(ismember(tracks(:,end-2),unique(labelMask(:))),:);
    
    if flag3d
        
        scat=scatter(handles.axes1,tracks(:,1),tracks(:,2),'rx');
    else
        scat=scatter(handles.axes1,tracks(:,1),tracks(:,2),'rx');
    end
    currentCentroids=tracks(:,[1,2,size(tracks,2)]);
    setappdata(handles.figure1,'currentCentroids',currentCentroids);
    
    hold(handles.axes1,'on')
    axis(handles.axes1,'equal');
    
    text(tracks(:,1),tracks(:,2),cellstr(num2str(tracks(:,end))),'VerticalAlignment'...
        ,'bottom', 'HorizontalAlignment','right','color',[1 1 1],'parent',handles.axes1);
  
    
    % plot boundaries
    B=bwboundaries(wormMask(:,:,zSlice));
    for i=1:length(B)
        b=B{i};
        plot(handles.axes1,b(:,2),b(:,1),'b')
    end
    hold(handles.axes1,'off')
    
else
    tracks=trackData(trackData(:,end-1)==iImage,:);
    scat3=getappdata(handles.figure1,'scat3');
    
    if ishandle(scat3)
        startP=repmat(mean(centroids),3,1);
        startP(:,2)=startP(:,2)+(-50:50:50)';
        
        %idx=kmeans(centroids(:,1:3),3,'start',startP);
        idx=3*ones(size(centroids,1),1);
        idx(centroids(:,2)<300)=2;
        
        idx(centroids(:,2)<200)=1;
        
        c=lines(max(idx));
        
        c=c(idx,:);
        %     c=c(idx==3);
        %     centroids=centroids(idx==3,:);
        offset=-mean(centroids(idx==2,1:2))+[150,250];
        centroids(idx==2,1:2)=bsxfun(@plus,centroids(idx==2,1:2),offset);
        
        set(scat3,'XData',centroids(:,2),'YData',centroids(:,1),...
            'ZData',centroids(:,3),'CData',c);
        %     hold(handles.axes1,'on');
        %         scatter3(handles.axes1,mean(centroids(idx==2,2)),...
        %             mean(centroids(idx==2,1)),mean(centroids(idx==2,3)));
        %          hold(handles.axes1,'off');
        
        
        ylim(handles.axes1,[0,size(wormMask,2)]);
        xlim(handles.axes1,[0,size(wormMask,1)])
        zlim(handles.axes1,[0,size(wormMask,3)*5])
    else
        idx=3*ones(size(centroids,1),1);
        idx(centroids(:,1)<300)=2;
        
        idx(centroids(:,1)<200)=1;
        
        c=lines(max(idx));
        c=c(idx,:);
        scat3=scatter3(handles.axes1,centroids(:,2),centroids(:,1),centroids(:,3),[],c);
        
        ylim(handles.axes1,[0,size(wormMask,2)]);
        xlim(handles.axes1,[0,size(wormMask,1)])
        zlim(handles.axes1,[0,size(wormMask,3)*5])
        
        setappdata(handles.figure1,'scat3',scat3);
    end
end

%  meanB=mean(cell2mat(B));
%
%     [meanBy,meanBx]=find(wormMask);
%     scatter(handles.axes1,mean(meanBx),mean(meanBy),'g')
%     yBase=sum(baseImg,2);
%     xBase=sum(baseImg,1);
%     yCM=dot(yBase,(1:length(yBase)))/sum(yBase);
%     xCM=dot(xBase,(1:length(xBase)))/sum(xBase);
%         scatter(handles.axes1,xCM,yCM,'gx')
%     hold(handles.axes1,'off')


if get(handles.showAll,'value')
    %show heatmap of all tracks
    activityMat=getappdata(handles.figure1,'activityMat');
    imagesc(activityMat,'parent',handles.axes2);
else
    
    %display point on axis 2
    displayIdx=get(handles.DisplayIdx,'data');
    plotIdx=[displayIdx{:,1}];
    plotIdx=plotIdx(~isnan(plotIdx) & plotIdx~=0);
    hold(handles.axes2,'off');
    
    switch get(handles.plotChannel,'value')
        case 1
            output=trackData(:,5);
        case 2
            output=trackData(:,4);
            
        case 3
            output=trackData(:,4)./trackData(:,5);
    end
    setappdata(handles.figure1,'output',output);
    
    switch get(handles.plotChannel2,'value')
        case 1
            output2= nan*ones(size(trackData(:,1)));
        case 2
            output2=trackData(:,5);
        case 3
            output2=trackData(:,4);
        case 4
            output2=trackData(:,4)./trackData(:,5);
        case 5
            output2=trackData(:,1);
    end
    setappdata(handles.figure1,'output',output);
    
    output(trackData(:,end-1)<startTime)=nan;
    %output=normalizeRange(output);
    %output=output/median(output);
    
    for i=1:length(plotIdx);
        idx=plotIdx(i);
        t=trackData((trackData(:,end)==idx),end-1);
        a=output((trackData(:,end)==idx));
        a2=output2((trackData(:,end)==idx));
        a=a(t>startTime);
        a2=a2(t>startTime);
        t=t(t>startTime);
        if normalizeFlag
            a=normalizeRange(smooth(a,smoothWindow))+i-1;
            a2=normalizeRange(smooth(a2,smoothWindow))+i-1;
            
        else
            
            a=(smooth(a,smoothWindow))+i-1;
            a2=(smooth(a2,smoothWindow))+i-1;
        end
        t=t*timeStep;
        plot(handles.axes2,t,a);
        hold(handles.axes2,'on');
        plot(handles.axes2,t,a2,'g');
    end
    
    
    
    
    subIdx=ismember(tracks(:,end),plotIdx);
    text(tracks(subIdx,1),tracks(subIdx,2),cellstr(num2str(tracks(subIdx,end))),'VerticalAlignment'...
        ,'bottom', 'HorizontalAlignment','right','color',[0 1 0],'parent',handles.axes1);
    
    
    hold(handles.axes2,'on');
    h=getappdata(0,'scatter');
    for iPlot=1:length(plotIdx);
        try
            delete(h(iPlot));
        catch
        end
        
        idx=plotIdx(iPlot);
        t=trackData((trackData(:,end)==idx),end-1);
        a=output((trackData(:,end)==idx));
        a=a(t>startTime);
        t=t(t>startTime);
        t=t*timeStep;
        if normalizeFlag
            a=normalizeRange(smooth(a,smoothWindow))+iPlot-1;
        else
            a=(smooth(a,smoothWindow))+iPlot-1;
        end
        a=a(t==(iImage*timeStep));
        
        if sum(a)
            a=a(1);
            h(iPlot)=scatter(handles.axes2,(iImage*timeStep),a,'r','fill');
        end
        
        hold(handles.axes2,'on');
    end
    setappdata(0,'scatter',h);
end
set(handles.currentFolder,'String',imFolder);


else
    if ~isempty(baseImg)
        ax1=imagesc(baseImg,'parent',handles.axes1);
        axis(handles.axes1,'equal')

    else
        scat3=scatter3(handles.axes1,centroids(:,2),centroids(:,1),centroids(:,3));
zlim(handles.axes1,[0,40])
ylim(handles.axes1,[0,600])
xlim(handles.axes1,[0,600])

    end
end
centerline=getappdata(handles.figure1,'centerline');
bfIdxLookup=getappdata(handles.figure1,'bfIdxLookup');

if ~isempty(centerline) && ~ isempty(bfIdxLookup)
    frameIdx=str2double(imName(6:10));
    bfIdx=bfIdxLookup(frameIdx);
    CLcurrent=[interp2(squeeze(centerline(:,1,:))',1:100,repmat(bfIdx,1,100))',...
        interp2(squeeze(centerline(:,2,:))',1:100,repmat(bfIdx,1,100))'];
    plot(handles.axes3,CLcurrent(:,2),CLcurrent(:,1));
    hold(handles.axes3,'on')
    scatter(handles.axes3,CLcurrent(1,2),CLcurrent(1,1),'xr')
    hold(handles.axes3,'off');
    xlim(handles.axes3,[1,1024]);
    ylim(handles.axes3,[1,1024]);
end
    


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
minDist=20;
params.mem=3;
params.good=100;
%params.excessive=3;
imFolder=getappdata(0,'imFolder');
rawImFolder=getappdata(0,'rawImFolder');

matFiles=dir([imFolder filesep '*.mat']);
setappdata(0,'matFiles',matFiles);
%imFiles=dir([imFolder filesep '*.tif']);
% setappdata(0,'imFiles',imFiles);

trackData=[];
trackIdx=0;
progressbar(0);
for imat=1:1:length(matFiles)
    trackIdx=trackIdx+1;
    x=load([imFolder filesep matFiles(imat).name]);
    if isfield(x,'centroids')
    centroids=x.centroids;
    Rintensities=x.Rintensities;
    Gintensities=x.Gintensities;
    Volume=x.Volume;
    time=x.time;
    wormMask=x.wormMask;
%    imFiles=x.imName;
%    baseImg= double(imread([rawImFolder{1} filesep imFiles]...
%       ,'tif')); %red
%         yBase=sum(baseImg,2);
%     xBase=sum(baseImg,1);
%     yCM=dot(yBase,(1:length(yBase)))/sum(yBase);
%     xCM=dot(xBase,(1:length(xBase)))/sum(xBase);
%         
%    centroids=bsxfun(@minus,centroids,(centroids'*Volume)'/sum(Volume));
    
    tracks=[centroids,Gintensities,Rintensities,trackIdx*ones(size(Gintensities))];
    trackData=[trackData;tracks];
    end
    progressbar((imat)/length(matFiles));
end

params.dim=size(centroids,2);
if any(trackData(:,3)<0)
trackData(:,3)=trackData(:,3)-min(trackData(:,3))+1;
end
%trackData=trackData(trackData(:,2)>200,:);
%trackData=trackData(trackData(:,end)<150,:);

for i=0:minDist*.5
    try
        display(['Min Dist is now ' num2str(minDist*.75^i)]);
        trackOutput=track(trackData,minDist*.75^i,params);
        break
    catch
        display(['reducing minDist by factor of ' num2str(i)]);
    end
 if i==minDist*.5
    display('Track Failed');
    trackOutput=ones(10);
    setappdata(handles.figure1,'trackOutput',trackOutput);
return

 end
        %    trackOutput=track(trackData,minDist*.75^i,params);

end

trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
%  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
[ trackIdx,ia,ib]=unique(trackOutput(:,end));
trackOutput(:,end)=ib;

nTracks=max(trackOutput(:,end));
nTime=max(trackOutput(:,end-1));
for iTrack=1:nTracks
    t=trackOutput((trackOutput(:,end)==iTrack),end-1);
    centroid=trackOutput((trackOutput(:,end)==iTrack),1:2);
    green=trackOutput((trackOutput(:,end)==iTrack),3);
    red=trackOutput((trackOutput(:,end)==iTrack),4);
    
    cellOutput(iTrack).time=t;
    cellOutput(iTrack).centroid=centroid;
    cellOutput(iTrack).green=green;
    cellOutput(iTrack).red=red;
    
    
end
setappdata(handles.figure1,'trackOutput',trackOutput);
save([imFolder filesep 'trackOutput'],'trackOutput','cellOutput')








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
trackData=getappdata(handles.figure1,'trackOutput');
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
