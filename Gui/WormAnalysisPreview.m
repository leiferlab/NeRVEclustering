function varargout = WormAnalysisPreview(varargin)
% WORMANALYSISPREVIEW MATLAB code for WormAnalysisPreview.fig
%      WORMANALYSISPREVIEW, by itself, creates a new WORMANALYSISPREVIEW or raises the existing
%      singleton*.
%
%      H = WORMANALYSISPREVIEW returns the handle to a new WORMANALYSISPREVIEW or the handle to
%      the existing singleton*.
%
%      WORMANALYSISPREVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORMANALYSISPREVIEW.M with the given input arguments.
%
%      WORMANALYSISPREVIEW('Property','Value',...) creates a new WORMANALYSISPREVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WormAnalysisPreview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WormAnalysisPreview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WormAnalysisPreview

% Last Modified by GUIDE v2.5 14-Apr-2017 14:29:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WormAnalysisPreview_OpeningFcn, ...
                   'gui_OutputFcn',  @WormAnalysisPreview_OutputFcn, ...
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


% --- Executes just before WormAnalysisPreview is made visible.
function WormAnalysisPreview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WormAnalysisPreview (see VARARGIN)

% Choose default command line output for WormAnalysisPreview
handles.output = hObject;

hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @plotter);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'SliderStep',[1,1]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WormAnalysisPreview wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function [hiResIdx,bfIdx]=getTimes(handles)
time=get(handles.slider1,'Value');
set(handles.currentTime,'String',num2str(round(time)));
bfAll=getappdata(handles.figure1,'bfAll');
hiResData=getappdata(handles.figure1,'hiResData');
hiResIdx=round(time);
bfTime=bfAll.frameTime;
bfIdx=interp1(bfTime,1:length(bfTime),hiResData.frameTime(hiResIdx),'nearest');
bfIdx(isnan(bfIdx))=1;

function  [redImage,greenImage,image0,image1]=getImages(handles)
alignments=getappdata(handles.figure1,'alignments');
rect1=alignments.S2AHiRes.rect1;
rect2=alignments.S2AHiRes.rect2;
dataFolder=getappdata(handles.figure1,'dataFolder');
[hiResIdx,bfIdx]=getTimes(handles);
Fid=getappdata(handles.figure1,'fileID');
sCMOSfile=dir([dataFolder filesep '*.dat']);
sCMOSfile=sCMOSfile.name;
sCMOSfile=[dataFolder filesep sCMOSfile ];
if isempty(Fid)
    Fid=fopen( sCMOSfile);
    setappdata(handles.figure1,'fileID',Fid);
elseif Fid<=0;
    Fid=fopen(sCMOSfile );
    setappdata(handles.figure1,'fileID',Fid);
end
[row,col]=getdatdimensions(sCMOSfile);
status=fseek(Fid,2*hiResIdx*row*col,-1);
temp=fread(Fid,row*col,'uint16',0,'l');
temp=(reshape(temp,row,col));
background=getappdata(0,'background');
temp=temp-background;
temp(temp<0)=0;

redImage=temp(rect1(2):rect1(4),rect1(1):rect1(3));
greenImage=temp(rect2(2):rect2(4),rect2(1):rect2(3));

cam0Obj=getappdata(handles.figure1,'cam0Obj');
cam1Obj=getappdata(handles.figure1,'cam1Obj');
if isempty(cam0Obj)
    image0=0;
    image1=0;
else
    image0=read(cam0Obj,bfIdx);
    image0=image0(:,:,1);
    image1=read(cam1Obj,bfIdx);
    image1=image1(:,:,1);
end

function CL=getCenterline(handles)
centerline=getappdata(handles.figure1,'centerline');
if ~isempty(centerline)
    [~,bfIdx]=getTimes(handles);
    CL=centerline(:,:,bfIdx);
else
    CL=[0 0; 0 0];
end




function plotter(hObject,eventdata)
handles=guidata(get(hObject,'Parent'));
[redImage,greenImage,image0,image1]=getImages(handles);
CL=getCenterline(handles);


            
if handles.viewMode.Value
    alignments=getappdata(handles.figure1,'alignments');
    Rsegment=alignments.S2AHiRes.Rsegment;
    
    switch handles.imageType.Value
        case 1
            bigImage=redImage;
        case 2
            tform=alignments.S2AHiRes.t_concord;
            bigImage=imwarp(greenImage,tform,'OutputView',Rsegment);
        case 3
    tform1=alignments.lowResFluor2BF.t_concord;
    tform2=alignments.Hi2LowResF.t_concord;
    tform1.T=tform1.T\tform2.T;
    bigImage=imwarp(image1,tform1,'OutputView',Rsegment);
        case 4
            tform=alignments.Hi2LowResF.t_concord;
            bigImage=imwarp(image0,tform,'OutputView',Rsegment);
    end
    
    bigPlot_handle=findobj(handles.bigPlot,'Type','Image');
    if ~isempty(bigPlot_handle)
        bigPlot_handle.CData=bigImage;
    else 
        imagesc(bigImage,'Parent',handles.bigPlot);
    end
    
    tform_lo2hi=alignments.lowResFluor2BF.t_concord;
    tform2=alignments.Hi2LowResF.t_concord;
    CL_hi=transformPointsInverse(tform_lo2hi,CL(:,[2,1]));
    CL_hi=transformPointsForward(tform2,CL_hi);
    plot_handle=findobj(handles.bigPlot,'Type','Line');
    if ~isempty(plot_handle)
        plot_handle.XData=CL_hi(:,1);
        plot_handle.YData=CL_hi(:,2);
    else
        hold(handles.bigPlot,'on');
        plot(handles.bigPlot,CL_hi(:,2),CL_hi(:,1),'r','LineWidth',4);
        hold(handles.bigPlot,'off');
    end
    
else
 
r40_handle=findobj(handles.r40,'Type','Image');
r40_handle.CData=redImage;

g40_handle=findobj(handles.g40,'Type','Image');
g40_handle.CData=greenImage;

dark10_handle=findobj(handles.dark10,'Type','Image');
dark10_handle.CData=image1;
plot_handle=findobj(handles.dark10,'Type','Line');
if ~isempty(plot_handle)
    plot_handle.XData=CL(:,2);
    plot_handle.YData=CL(:,1);
    
else
    hold(handles.dark10,'on');
    plot(handles.dark10,CL(:,2),CL(:,1),'r','LineWidth',2);
    hold(handles.dark10,'off');
end


fluor10_handle=findobj(handles.fluor10,'Type','Image');
fluor10_handle.CData=image0;
end

% --- Outputs from this function are returned to the command line.
function varargout = WormAnalysisPreview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%select image folder
display('Select a data folder');
display(['For an example working folder, select the folder at '...
    '/tigress/LEIFER/PanNeuronal/testing_sets/BrainScanner20161031_111303'])
mostRecent=getappdata(0,'mostRecent');
if isempty(mostRecent)
    dataFolder=uipickfiles('Prompt', 'Select the data Folder (BrainScanner)' );
else
    dataFolder=uipickfiles('filterspec', mostRecent,...
        'Prompt', 'Select the data Folder (BrainScanner)');
end
dataFolder=dataFolder{1};
mostRecent=fileparts(dataFolder);


setappdata(0,'mostRecent',mostRecent);
set(handles.currentFolder,'String',dataFolder);

%load avi movies
lowMagFolder=dir([dataFolder filesep 'LowMag*']);
lowMagFolder=[dataFolder filesep lowMagFolder.name];
cam0File=[lowMagFolder filesep 'cam0.avi'];
cam1File=[lowMagFolder filesep 'cam1.avi'];
try
    
    cam0Obj=VideoReader(cam0File);
    cam1Obj=VideoReader(cam1File);
    setappdata(handles.figure1,'cam0Obj',cam0Obj);
    setappdata(handles.figure1,'cam1Obj',cam1Obj);
    
    set(handles.lowMagStatus,'String', 'LowMag Present');
    handles.lowMagStatus.BackgroundColor=[0,.94,0];
    set(handles.lowMagStatus,'Enable','off');
catch me
    set(handles.lowMagStatus,'Enable','on');
    set(handles.lowMagStatus,'String', 'LowMag Missing');
    handles.lowMagStatus.BackgroundColor=[.94,0,0];
    setappdata(handles.figure1,'cam0Obj',[]);
    setappdata(handles.figure1,'cam1Obj',[]);
end


%load timing data
try
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder);
        set(handles.timingStatus,'String', 'Timing Good');
        handles.timingStatus.BackgroundColor=[0,.94,0];
    set(handles.timingStatus,'Enable','off');
catch me
        set(handles.timingStatus,'String', me.identifier);
        handles.timingStatus.BackgroundColor=[.94,0,0];
        handles.timingStatus.Callback={@popupForCallbacks,me.message};
        set( handles.timingStatus,'Enable','on');
end


%load alignments
try
    alignments=load([dataFolder filesep 'alignments.mat']);
    alignments=alignments.alignments;
    setappdata(handles.figure1,'alignments',alignments);
    
    if isfield(alignments,'background')
        background=alignments.background;
        set(handles.alignmentStatus,'String', 'Alignments Present');
        handles.alignmentStatus.BackgroundColor=[0,.94,0];
    else
        background=0;
        set(handles.centerlineStatus,'String', 'Background Missing');
        handles.alignmentStatus.BackgroundColor=[0,.94,0];
    end
    set(handles.alignmentStatus,'Enable','off');
catch me
    set(handles.centerlineStatus,'String', 'Alignments Missing');
    handles.alignmentStatus.BackgroundColor=[0,.94,0];
    set(handles.alignmentStatus,'Enable','on');
    handles.alignmentStatus.Callback={@popupForCallbacks,me.message};
end

setappdata(handles.figure1,'background',background);



%load centerlines
try
    [centerline,~, ~,~]=loadCLBehavior(dataFolder);
    setappdata(handles.figure1,'centerline',centerline);
    set(handles.centerlineStatus,'String', 'Centerline Present');
    handles.centerlineStatus.BackgroundColor=[0,.94,0];
    set(handles.centerlineStatus,'Enable','off');
catch me
    set(handles.centerlineStatus,'String', 'Centerline Missing');
    handles.centerlineStatus.BackgroundColor=[.94,0,0];
    setappdata(handles.figure1,'centerline',[]);
    set(handles.centerlineStatus,'Enable','on');
    handles.centerlineStatus.Callback={@popupForCallbacks,me.message};

end


if exist('me','var')
    return
end
setappdata(handles.figure1,'hiResData',hiResData);
setappdata(handles.figure1,'bfAll',bfAll);
setappdata(handles.figure1,'dataFolder',dataFolder);

%setting slider parameters
set(handles.slider1,'Min',1)
set(handles.slider1,'Value',1)
setappdata(handles.figure1,'cursorTarget', 1);

maxFrame=length(hiResData.stackIdx);
minFrame=1;
set(handles.maxTime,'String',num2str(maxFrame));
set(handles.minTime,'String',num2str(minFrame));
set(handles.slider1,'max',maxFrame);
setappdata(handles.figure1,'currentFrame',1);
set(handles.slider1,'value',1);

[redImage,greenImage,image0,image1]=getImages(handles);
imagesc(redImage,'Parent',handles.r40);
imagesc(greenImage,'Parent',handles.g40);
imagesc(image1,'Parent',handles.dark10);
imagesc(image0,'Parent',handles.fluor10);
axis(handles.r40,'off');
axis(handles.g40,'off');
axis(handles.dark10,'off');
axis(handles.fluor10,'off');
axis(handles.bigPlot,'off');



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


% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value=handles.slider1.Value;
handles.slider1.Value=max(1,value-1);
plotter(hObject,handles);


function currentTime_Callback(hObject, eventdata, handles)
% hObject    handle to currentTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentTime as text
%        str2double(get(hObject,'String')) returns contents of currentTime as a double
value=str2double(get(hObject,'String'));
set(handles.slider1,'Value',value);
plotter(handles.slider1,eventdata);


% --- Executes during object creation, after setting all properties.
function currentTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value=handles.slider1.Value;
handles.slider1.Value=min(handles.slider1.Max,value+1);
plotter(hObject,handles);


% --- Executes on button press in go2flash.
function go2flash_Callback(hObject, eventdata, handles)
% hObject    handle to go2flash (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hiResData=getappdata(handles.figure1,'hiResData');
flashLoc=hiResData.flashLoc;
currentTime=get(handles.slider1,'Value');
nextTime=flashLoc(find(flashLoc>currentTime,1,'first'));
if any(nextTime)
    handles.currentTime.String=num2str(nextTime);
    currentTime_Callback(handles.currentTime, eventdata, handles);
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


% --- Executes on button press in viewMode.
function viewMode_Callback(hObject, eventdata, handles)
% hObject    handle to viewMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of viewMode
if get(hObject,'Value')
    handles.bigPlot.Visible= 'on';
    hObject.BackgroundColor=[.94,0,0];
    uistack(handles.bigPlot,'top')
    axis(handles.bigPlot,'off');
    set(handles.imageType,'Visible','on')

else
    handles.bigPlot.Visible= 'off';
    delete(handles.bigPlot.Children);
    hObject.BackgroundColor=[.94,.94,.94];
    uistack(handles.bigPlot,'bottom')
    set(handles.imageType,'Visible','off')
end


plotter(hObject,handles);


% --- Executes on selection change in imageType.
function imageType_Callback(hObject, eventdata, handles)
% hObject    handle to imageType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imageType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageType
plotter(hObject,handles);

% --- Executes during object creation, after setting all properties.
function imageType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, evnt, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmp(evnt.Key,'rightarrow')|| strcmp(evnt.Key,'d')
    right_Callback(handles.slider1,evnt,handles);
    %Backward
elseif strcmp(evnt.Key,'backspace') || strcmp(evnt.Key,'leftarrow')|| strcmp(evnt.Key,'a')
    left_Callback(handles.slider1,evnt,handles);
elseif strcmp(evnt.Key,'shift');
    image_val=handles.imageType.Value;
    image_val=image_val+1;
    if image_val>4, image_val=1; end
    handles.imageType.Value=image_val;
    plotter(handles.slider1,evnt);
end


% --- Executes on button press in instructions.
function instructions_Callback(hObject, eventdata, handles)
% hObject    handle to instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function popupForCallbacks(hObject,eventdata,message)

%make a pop up message box, for showing information about an error usually
msgbox(message,'Error!','error')
