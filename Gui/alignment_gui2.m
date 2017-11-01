function varargout = alignment_gui(varargin)
% ALIGNMENT_GUI MATLAB code for alignment_gui.fig
% This gui is made to help streamline the process of creating alignment
% files for the Leifer lab whole brain imaging set up. The setup produces 4
% video feeds from 3 cameras to track and measure neural activity from
% moving C. elegans. It is therefore necessary to spatially align all of
% the video feeds. This is done by taking several images of different
% fields of view of tetraspek beads. The bead images are taken with the
% same whole brain imaging programs used to image the worm. This program
% will take those recordings, sepearte each individual field of view, and
% allow the user to manually click on the bead in each field of view for
% each of the 4 images. The clicked coordinates can then be used to create
% an alignment.mat file used in the wholebrain analysis pipeline. 


%      ALIGNMENT_GUI, by itself, creates a new ALIGNMENT_GUI or raises the existing
%      singleton*.
%
%      H = ALIGNMENT_GUI returns the handle to a new ALIGNMENT_GUI or the handle to
%      the existing singleton*.
%
%      ALIGNMENT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGNMENT_GUI.M with the given input arguments.
%
%      ALIGNMENT_GUI('Property','Value',...) creates a new ALIGNMENT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alignment_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alignment_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alignment_gui

% Last Modified by GUIDE v2.5 01-Nov-2017 10:24:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @alignment_gui_OpeningFcn, ...
    'gui_OutputFcn',  @alignment_gui_OutputFcn, ...
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


% --- Executes just before alignment_gui is made visible.
function alignment_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alignment_gui (see VARARGIN)

% Choose default command line output for alignment_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
%set(src, 'Pointer', 'crosshair'); % Optional
pan off % Panning will interfere with this code

% UIWAIT makes alignment_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alignment_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
display([
    'For testing purposes, use della VNC and select one of the folders in'...
    '/tigress/LEIFER/PanNeuronal/testing_sets/alignments_test'])
mostRecent=getappdata(0,'mostRecent');
imFolder=uipickfiles('filterspec',mostRecent,...
    'Prompt','Select folder containing alignment .avi and .dat files');
imFolder=imFolder{1};
setappdata(0,'mostRecent',imFolder);
setappdata(handles.figure1,'imFolder',imFolder);

[dat_images,avi0_images,avi1_images]=ExtractAlignmentImagesFromVideo(imFolder);


red_images=dat_images(1:512,:,:);
green_images=dat_images(513:end,:,:);
nFrames=min(size(avi0_images,3),size(red_images,3));
all_images=repmat(struct(),1,nFrames);
setappdata(handles.figure1,'nFrames',nFrames);
for iFrame=1:nFrames
    all_images(iFrame).red=red_images(:,:,iFrame);
    all_images(iFrame).green=green_images(:,:,iFrame);
    all_images(iFrame).fluor10=avi0_images(:,:,iFrame);
    all_images(iFrame).dark10=avi1_images(:,:,iFrame);
end
setappdata(handles.figure1,'all_images',all_images)

axis_list={handles.red_image,...
    handles.green_image,...
    handles.fluor10,...
    handles.dark10...
    };

setappdata(handles.figure1,'axisList',axis_list);
setappdata(handles.figure1,'zoom_flag',0);

showplots(handles)


% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFrame=str2num(get(handles.currentFrame,'String'));
newFrame=max(currentFrame-1,1);
set(handles.currentFrame,'String',num2str(newFrame));
showplots(handles)

% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFrame=str2num(get(handles.currentFrame,'String'));
maxFrame=getappdata(handles.figure1,'nFrames');
newFrame=min(currentFrame+1,maxFrame);
set(handles.currentFrame,'String',num2str(newFrame));
showplots(handles)

% --- Display 4 plots and annotated points
function showplots(handles)
% handles    structure with handles and user data (see GUIDATA)
all_images=getappdata(handles.figure1,'all_images');
currentFrame=str2num(get(handles.currentFrame,'String'));
images=all_images(currentFrame);
axis_list=getappdata(handles.figure1,'axisList');

image_names=fieldnames(images);

demoState=[...
    get(handles.redView,'Value'),...
    get(handles.greenView,'Value'),...
    get(handles.fluorView,'Value'),...
    get(handles.darkView,'Value')];
label_list={handles.nRed,handles.nGreen,handles.nFluor,handles.nDark};
for iAxes=1:4
    
    current_axis=axis_list{iAxes};
    im_handles=findobj(current_axis,'Type','Image');
    if isempty(im_handles)
        imagesc(images.(image_names{iAxes}),'parent',current_axis);
        axis(current_axis,'off')
    else
        im_handles.CData=images.(image_names{iAxes});
    end
    delete(findobj(current_axis,'Type','Scatter'))
    delete(findobj(current_axis,'Type','Text'))
    
    
    hold(current_axis,'on')
    points=getappdata(current_axis,'points');
    if any(points)
        scatter(current_axis,points(:,1),points(:,2),'xr')
        text(points(:,1),...
            points(:,2),...
            cellstr(num2str([1:size(points,1)]')),...
            'parent',current_axis,...
            'color','white',...
            'HorizontalAlignment','right');
    end

    
    if any(demoState) && find(demoState)~=iAxes
        demoPoints=getappdata(current_axis,'demoPoints');
        scatter(current_axis,demoPoints(:,1),demoPoints(:,2),'og')
        text(demoPoints(:,1),...
            demoPoints(:,2),...
            cellstr(num2str([1:size(demoPoints,1)]')),...
            'parent',current_axis,...
            'color','magenta',...
            'HorizontalAlignment','left');
    end
    
    hold(current_axis,'off')
    set(label_list{iAxes},'String',num2str((size(points,1))));
    
end





% --- Executes on button press in saveTransformations.
function saveTransformations_Callback(hObject, eventdata, handles)
% hObject    handle to saveTransformations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get original folder and make output names
imFolder=getappdata(handles.figure1,'imFolder');
[~,fileName]=fileparts(imFolder);
HiResS2A_name=[fileName 'HiResS2A'];
LowResB2F_name=[fileName 'LowResB2F'];
HiResS2LoResF_name=[fileName 'HiResS2LoResF'];

[~,hostname]=system('hostname');
if strfind(hostname,'tigressdata')
    alignmentFolder='/tigress/LEIFER/communalCode/3dbrain/registration/';
    backgroundLocation='/tigress/LEIFER/communalCode/3dbrain/Background';
else
alignmentFolder='Y:\CommunalCode\3dbrain\registration\';
backgroundLocation='Y:\CommunalCode\3dbrain\background\';
end

alignments=makeTransformation(handles);
setappdata(handles.figure1,'alignments',alignments)
% if there's a background image, load it as well into alignments.
display('select a camera background image for this size himag video');

backgroundImage=uipickfiles('FilterSpec',backgroundLocation,...
    'Prompt','Select Background HiMag Image File');

if iscell(backgroundImage)
    backgroundImage=load(backgroundImage{1});
    backgroundImage=backgroundImage.backgroundImage;
else
    backgroundImage=0;
end
alignments.background=backgroundImage;

% save 
S2AHiRes=alignments.S2AHiRes;
Hi2LowResF=alignments.Hi2LowResF;
lowResFluor2BF=alignments.lowResFluor2BF;

save([alignmentFolder HiResS2A_name],'-struct'...
    ,'S2AHiRes')
save([alignmentFolder HiResS2LoResF_name],'-struct'...
    ,'Hi2LowResF')
save([alignmentFolder LowResB2F_name],'-struct'...
    ,'lowResFluor2BF')
save([imFolder filesep 'alignments'],'alignments');
display('Saved! Copy the alignment.mat folder into all the BrainScanner folders.')


function alignments=makeTransformation(handles)

axis_list=getappdata(handles.figure1,'axisList');

%get all the coordinates
points_cell=cell(1,4);
for iAxes=1:4
    current_axis=axis_list{iAxes};
    points=getappdata(current_axis,'points');
    points_cell{iAxes}=points;
end

%get images for size purposes
all_images=getappdata(handles.figure1,'all_images');
images=all_images(1);

% save in a consistent way, will change eventually to clean this up

%Select channel for segmentation depending on indicator (from menu in gui)
%1: not selected, defaulting to green
%2: green
%3: red
if handles.calciumchannel==1 || handles.calciumchannel==2
    S2AHiRes.rect1=[1,1,512,512];
    S2AHiRes.rect2=[1,513,512,512];
elseif handles.calciumchannel==3
    S2AHiRes.rect2=[1,1,512,512];
    S2AHiRes.rect1=[1,513,512,512];
end

S2AHiRes.Aall=points_cell{2};
S2AHiRes.Sall=points_cell{1};
S2AHiRes.t_concord = ...
    fitgeotrans(S2AHiRes.Aall,S2AHiRes.Sall,'projective');
S2AHiRes.Rsegment = imref2d(size(images.red));

Hi2LowResF.Aall=points_cell{3};
Hi2LowResF.Sall=points_cell{1};
Hi2LowResF.t_concord = ...
    fitgeotrans(Hi2LowResF.Aall,Hi2LowResF.Sall,'projective');
Hi2LowResF.Rsegment = imref2d(size(images.red));

lowResFluor2BF.Aall=points_cell{3};
lowResFluor2BF.Sall=points_cell{4};
lowResFluor2BF.t_concord = ...
    fitgeotrans(lowResFluor2BF.Aall,lowResFluor2BF.Sall,'projective');
lowResFluor2BF.Rsegment = imref2d(size(images.dark10));

% if you select them individually, bundle them and save it for later use
alignments.lowResFluor2BF=lowResFluor2BF;
alignments.S2AHiRes=S2AHiRes;
alignments.Hi2LowResF=Hi2LowResF;






% --- Executes on button press in loadAlignment.
function loadAlignment_Callback(hObject, eventdata, handles)
% hObject    handle to loadAlignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imFolder=getappdata(handles.figure1,'imFolder');
alignment_file=uipickfiles('filterspec',imFolder,...
    'prompt', 'Select Alignment File');

alignment_file=alignment_file{1};
alignments=load(alignment_file);
alignments=alignments.alignments;
setappdata(handles.figure1,'alignment',alignments);

axis_list=getappdata(handles.figure1,'axisList');

%get all the coordinates, adding points from old alignments top whats there
%already. 
points_cell=cell(1,4);
for iAxes=1:4
    current_axis=axis_list{iAxes};
    points=getappdata(current_axis,'points');
    switch iAxes
        case 1
            newPoints=alignments.S2AHiRes.Sall;
        case 2
            newPoints=alignments.S2AHiRes.Aall;
        case 3
            newPoints=alignments.Hi2LowResF.Aall;
        case 4
            newPoints=alignments.lowResFluor2BF.Sall;
    end
    points=[newPoints;points];
    setappdata(current_axis,'points',points);

end


% --- Executes on mouse press over axes background.
function red_image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to red_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on mouse press over axes background.
function green_image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to green_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function fluor10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fluor10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function dark10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dark10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Execute on mouse press on one of the axes. Left click adds new point
% at cursor location, right click removes most recent point from that axis.
% also zooms if the zoomflag is on
function getMousePositionOnImage(src, event)

handles = guidata(src);
axis_list=getappdata(handles.figure1,'axisList');
for iAxes=1:4
    current_axis=axis_list{iAxes};
    cursorPoint = get(current_axis, 'CurrentPoint');
    curX = cursorPoint(1,1);
    curY = cursorPoint(1,2);
    
    xLimits = get(current_axis, 'xlim');
    yLimits = get(current_axis, 'ylim');
    
    if (curX > min(xLimits) && curX < max(xLimits) &&...
            curY > min(yLimits) && curY < max(yLimits))
         zoom_flag=getappdata(handles.figure1,'zoom_flag');
        if zoom_flag
            if strcmp(src.SelectionType,'normal')
                zoomcenter(current_axis,curX,curY,2);
            else
                zoomcenter(current_axis,curX,curY,.5);
            end
        else
            if strcmp(src.SelectionType,'normal')
                addPoint(current_axis,[curX,curY]);
            else
                removePoint(current_axis);
            end
            break
        end
    end
end


showplots(handles)

function addPoint(current_axis,point)
old_points=getappdata(current_axis,'points');
setappdata(current_axis,'points',[old_points;point]);

function removePoint(current_axis)
old_points=getappdata(current_axis,'points');
setappdata(current_axis,'points',old_points(1:end-1,:));


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key
    case 'a'
        left_Callback(handles.left, eventdata, handles)
    case 'd'
        right_Callback(handles.right, eventdata, handles)
    case 'z'
        handles.figure1.Pointer='circle';
        setappdata(handles.figure1,'zoom_flag',1);
end




% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key

    case 'z'
        handles.figure1.Pointer='arrow';
        setappdata(handles.figure1,'zoom_flag',0);
end

function zoomcenter(varargin)
%ZOOMCENTER Zoom in and out of a specifeid point on a 2-D plot.
% ZOOMCENTER(X,Y) zooms the current axis on the point (X,Y) by a factor of 2.5.
% ZOOMCENTER(X,Y,FACTOR) zooms the current axis on the point (X,Y) by FACTOR.
%
% ZOOMCENTER(AX,...) zooms on the specified axis
%
% Example:
% line
% zoomcenter(.5, .5, 10)
%
% line
% zoomcenter(.7, .3, .5)

nin = nargin;
if ishandle(varargin{1})
 ax = varargin{1};
 varargin = varargin(2:end);
 nin = nin-1;
else
 ax = gca;
end
if nin<2
 error('ZOOMCENTER requires specifying both X and Y');
else
 x = varargin{1};
 y = varargin{2};
end
if nin==3
 factor = varargin{3};
else
 factor = 2.5;
end
im_handles=findobj(ax,'Type','Image');

imsize=size(im_handles.CData);

cax = axis(ax);
daxX = (cax(2)-cax(1))/factor(1)/2;
daxY = (cax(4)-cax(3))/factor(end)/2;
cax_new=[x+[-1 1]*daxX y+[-1 1]*daxY];
cax_new(cax_new<1)=1;
cax_new(2)=min(cax_new(2),imsize(1));
cax_new(4)=min(cax_new(4),imsize(2));

xlim(ax,cax_new([1,2]));
ylim(ax,cax_new([3,4]));


% --- Executes on button press in redView.
function redView_Callback(hObject, eventdata, handles)
% hObject    handle to redView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of redView
if get(hObject,'Value')
    setDemoView(hObject,handles,1);
    hObject.BackgroundColor=[1 .4 .4];
    handles.greenView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.darkView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.fluorView.BackgroundColor=[0.9400 0.9400 0.9400];
else
    hObject.BackgroundColor=[0.9400 0.9400 0.9400];
end
showplots(handles)


% --- Executes on button press in fluorView.
function fluorView_Callback(hObject, eventdata, handles)
% hObject    handle to fluorView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fluorView
if get(hObject,'Value')
    setDemoView(hObject,handles,3);
    hObject.BackgroundColor=[1 .4 .4];
    handles.redView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.darkView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.greenView.BackgroundColor=[0.9400 0.9400 0.9400];
else
    hObject.BackgroundColor=[0.9400 0.9400 0.9400];
end
showplots(handles)


% --- Executes on button press in greenView.
function greenView_Callback(hObject, eventdata, handles)
% hObject    handle to greenView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of greenView
if get(hObject,'Value')
    setDemoView(hObject,handles,2);
    hObject.BackgroundColor=[1 .4 .4];
    handles.redView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.darkView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.fluorView.BackgroundColor=[0.9400 0.9400 0.9400];
else
    hObject.BackgroundColor=[0.9400 0.9400 0.9400];
end
showplots(handles)



% --- Executes on button press in darkView.
function darkView_Callback(hObject, eventdata, handles)
% hObject    handle to darkView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of darkView
if get(hObject,'Value')
    setDemoView(hObject,handles,4);
    hObject.BackgroundColor=[1 .4 .4];
    handles.redView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.greenView.BackgroundColor=[0.9400 0.9400 0.9400];
    handles.fluorView.BackgroundColor=[0.9400 0.9400 0.9400];
    
else
    hObject.BackgroundColor=[0.9400 0.9400 0.9400];
end
showplots(handles)



function setDemoView(hObject,handles,axisID)

%reset all other switches
set(handles.redView,'Value',0);
set(handles.greenView,'Value',0);
set(handles.fluorView,'Value',0);
set(handles.darkView,'Value',0);
set(hObject,'Value',~get(hObject,'Value'));

axis_list=getappdata(handles.figure1,'axisList');

%get all the coordinates
points_cell=cell(1,4);
for iAxes=1:4
    current_axis=axis_list{iAxes};
    points=getappdata(current_axis,'points');
    points_cell{iAxes}=points;
end

for iAxes=1:4
    current_axis=axis_list{iAxes};
    t=fitgeotrans( points_cell{axisID}, points_cell{iAxes},'projective');
    demoPoints=transformPointsForward(t, points_cell{axisID});
    setappdata(current_axis,'demoPoints',demoPoints)
end


% --- Executes during object creation, after setting all properties.
function calciumchannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calciumchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_calciumChannel(hObject,eventdata,handles)
    %S2AHiRes.rect1=[1,1,512,512];
    %S2AHiRes.rect2=[1,513,512,512];
    vars = guidata(hObject);
    handles.calciumchannel = vars.calciumchannel.Value;
    %handles.rect1=
