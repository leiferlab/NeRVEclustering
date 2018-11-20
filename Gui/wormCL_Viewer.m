function varargout = wormCL_Viewer(varargin)
% WORMCL_VIEWER MATLAB code for wormCL_Viewer.fig
%      WORMCL_VIEWER, by itself, creates a new WORMCL_VIEWER or raises the existing
%      singleton*.
%
%      H = WORMCL_VIEWER returns the handle to a new WORMCL_VIEWER or the handle to
%      the existing singleton*.
%
%      WORMCL_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORMCL_VIEWER.M with the given input arguments.
%
%      WORMCL_VIEWER('Property','Value',...) creates a new WORMCL_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wormCL_Viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wormCL_Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wormCL_Viewer

% Last Modified by GUIDE v2.5 01-Mar-2018 16:05:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wormCL_Viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @wormCL_Viewer_OutputFcn, ...
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


% --- Executes just before wormCL_Viewer is made visible.
function wormCL_Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wormCL_Viewer (see VARARGIN)

% Choose default command line output for wormCL_Viewer
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @showImage);


setappdata(handles.figure1,'holdaxes',false);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'Max',2000);
%set defaults 
set(handles.colorMap,'Value',3)
set(handles.transpose,'Value', true)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wormCL_Viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wormCL_Viewer_OutputFcn(hObject, eventdata, handles) 
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
setappdata(0,'mostRecent',fileparts(currentData));
set(handles.currentFolder,'String',currentData);

Fid= VideoReader(currentData);
nFrames=Fid.NumberOfFrames;
setappdata(handles.figure1,'Fid',Fid);


setappdata(handles.figure1,'maxC',0);
setappdata(handles.figure1,'nFrames',nFrames);
set(handles.slider1,'Value',1);
set(handles.slider1,'Max',nFrames);
set(handles.maxSlider,'String',num2str(nFrames));
set(handles.minSlider,'String','1');
setappdata(handles.figure1,'MaxValue',nFrames);


%Also get CLworkspace if its present, useful for subtracting backgrounds
parentFolder=fileparts(currentData);
CLworkspace_file=[parentFolder filesep 'CLworkspace.mat'];
if exist(CLworkspace_file,'file')
    display('Workspace found, loading workspace')
    CLworkspace=load(CLworkspace_file);
    try 
        background=CLworkspace.mean_bf_all;
        frame_bg_lvl=CLworkspace.frame_bg_lvl;
        set(handles.backgroundSubtract,'enable','on');
        cline_para=CLworkspace.cline_para;
        setappdata(handles.figure1,'cline_para',cline_para);
    catch
        background=[];
        frame_bg_lvl=[];
        set(handles.backgroundSubtract,'enable','off');
    end
else
    background=[];
    frame_bg_lvl=[];
    set(handles.backgroundSubtract,'enable','off');
end
    display('Load tip file if present, otherwise, cancel')
    tip_file=[parentFolder filesep 'tip_coodinates.mat'];
    if exist(tip_file,'file')
        tips=load(tip_file);
        display(' Tip file found, loading tips');
        % update tip variable in CLWorkspace
        CLworkspace.tips = tips
    try
    tips=processTips(CLworkspace);
    setappdata(handles.figure1,'tips',tips);
    catch
    end
    
    
        
        % save tips for GUI
        setappdata(handles.figure1,'tips',tips);
    else
        tips=[];
        display('No tips found!');
    end
setappdata(handles.figure1,'background',background);
setappdata(handles.figure1,'frame_bg_lvl',frame_bg_lvl);


showImage(hObject)


function showImage(hObject,eventdata)

handles=guidata(get(hObject,'Parent'));
frameNumber=get(handles.slider1,'Value');
frameNumber=max(1,round(frameNumber));
Fid=getappdata(handles.figure1,'Fid');
CL_set=getappdata(handles.figure1,'centerline');
CLoffset=getappdata(handles.figure1,'CLoffset');
CLnumber=frameNumber+CLoffset;
    C=read(Fid,frameNumber);
    C=C(:,:,1);


if handles.backgroundSubtract.Value
    C=processImage(handles);
end

maxC=getappdata(handles.figure1,'maxC');
setappdata(handles.figure1,'maxC',max(max(C(:)),maxC));

h=imagesc(C,'Parent',handles.axes1);
set(handles.currentFrame,'String',num2str(frameNumber));
switch get(handles.colorMap,'Value')
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
        
    otherwise
end
hold(handles.axes1,'on')
colors='rbgcy';
for iCL=1:length(CL_set)
    centerline=CL_set{iCL};
if  get(handles.transpose,'Value')
    X=centerline(:,2,CLnumber);
    Y=centerline(:,1,CLnumber);
else
    X=centerline(:,1,CLnumber);
    Y=centerline(:,2,CLnumber);
end
if CLnumber<size(centerline,3)
plot(handles.axes1,X, Y,colors(iCL));
end
plot(handles.axes1,X(end,:), Y(end,1),'og');
plot(handles.axes1,X(1,:), Y(1,1),'ob');
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
showImage(hObject)


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
showImage(hObject)
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


if isempty(ender);;
    ender=1;
else
    ender=ender+1;
end

setappdata(handles.figure1,'ender',ender);


tiffwrite(imageName,single(currentImage),'tif',0);

currentFrame=get(handles.slider1,'Value');
set(handles.slider1,'Value',min(currentFrame+1,get(handles.slider1,'max')));
showImage(hObject)

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
for i_file=1:length(CLfile)
CLfiletemp=CLfile{i_file};
centerline=load(CLfiletemp);
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
CL_set{i_file}=centerline;
end

setappdata(handles.figure1,'centerline',CL_set);


% --- Executes on button press in transpose.
function transpose_Callback(hObject, eventdata, handles)
% hObject    handle to transpose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transpose


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imcontrast(handles.axes1)


% --- Executes on button press in loadTips.
function loadTips_Callback(hObject, eventdata, handles)
% hObject    handle to loadTips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFolder=get(handles.currentFolder,'String');
currentFolder=fileparts(currentFolder);

tip_file=uipickfiles('Filterspec', currentFolder);
tip_file=tip_file{1};
tip_data=load(tip_file);
head_pts=tip_data.head_pts;
tail_pts=tip_data.tail_pts;

setappdata(handles.figure1,'head_pts',head_pts);
setappdata(handles.figure1,'tail_pts',tail_pts);



% --- Executes on button press in fixTips.
function fixTips_Callback(hObject, eventdata, handles)
% hObject    handle to fixTips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in initCL.
function initCL_Callback(hObject, eventdata, handles)
% hObject    handle to initCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentFrame=get(handles.slider1,'Value');
currentFrame = max(1,round(currentFrame));

[y,x] = getpts(handles.axes1);
P = [x(:) y(:)];

O(:,1)=interp([P(:,1)],10);
O(:,2)=interp([P(:,2)],10);
nPoints = 100;

dis=[0;cumsum(sqrt(sum((O(2:end,:)-O(1:end-1,:)).^2,2)))];

K(:,1) = interp1(dis,O(:,1),linspace(0,dis(end),nPoints));
K(:,2) = interp1(dis,O(:,2),linspace(0,dis(end),nPoints));
centerline=getappdata(handles.figure1,'centerline');
centerline=centerline{1};
centerline(:,:,currentFrame)=K;
setappdata(handles.figure1,'centerline',{centerline});
 showImage(handles.slider1,eventdata)


% --- Executes on button press in autoCL.
function autoCL_Callback(hObject, eventdata, handles)
% hObject    handle to autoCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%fit a certain number of frames with autoCL
    %get centerline
    for i=1:str2double(get(handles.stepSize,'String'))
    cline_para=getappdata(handles.figure1,'cline_para');

    frameNumber=get(handles.slider1,'Value');
    frameNumber=max(1,round(frameNumber));
    CL_set=getappdata(handles.figure1,'centerline');
    if i==1
    cl_old=CL_set{1}(:,[2,1],frameNumber);
    else
        cl_old=cl;
    end
    tips=getappdata(handles.figure1,'tips');
        if ~isempty(tips)
        cline_para.head_pt=tips.head_pts(frameNumber,:);
        cline_para.tail_pt=tips.tail_pts(frameNumber,:);

        cline_para.stretching_force_factor=[0.05 0.05];
        cline_para.stretch_ends_flag=1;
        % Monika added this
        %cline_para.CLbeta=20;
        %cline_para.gamma=60;
        %cline_para.heat=10;
        %cline_para.refL=20.5;
        cline_para.gradient_force=1;
        cline_para.endRepulsion=.01;
        %end Monika

    end
    cline_para.stretching_force_factor=[0.1 0.1];

    cline_para.memForce=0.1;
    C=processImage(handles);
    cline_para.CLalpha=1;
    %cline_para.CLbeta = 1;
    cline_para.gradient_force=1;
    cline_para.endRepulsion=.01;
    cline_para.refSpring=0.00;
    cline_para.showFlag=0;
    % Monika
    %cline_para.endKappa=0.05;
    
    tip_image=C;
    cl_old=distanceInterp(cl_old(1:end,:),100);
    [cl,Is,Eout]=ActiveContourFit_wormRef4(...
        C,tip_image, cline_para, cl_old,1,[0 0]);
     showImage(handles.slider1,eventdata)

    CL_set{1}(:,[2,1],frameNumber)=cl;
    setappdata(handles.figure1,'centerline',CL_set);
    hold(handles.axes1,'on')
    plot(handles.axes1,cl(:,1),cl(:,2))
    hold(handles.axes1,'off');
   
    handles.slider1.Value=handles.slider1.Value+1;
   

end

% --- Executes on button press in backgroundSubtract.
function backgroundSubtract_Callback(hObject, eventdata, handles)
% hObject    handle to backgroundSubtract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showImage(handles.slider1,eventdata)

% Hint: get(hObject,'Value') returns toggle state of backgroundSubtract

function tips=processTips(CLworkspace)
    %if the tips are present process them and output the head and tail pts
    
    head_pts=CLworkspace.tips.head_pts;
    tail_pts=CLworkspace.tips.tail_pts;
    
    tip_mat_length=min(size(head_pts,1),size(tail_pts,1));
    head_pts=head_pts(1:tip_mat_length,:);
    tail_pts=tail_pts(1:tip_mat_length,:);
    
    same_pt=all(head_pts==tail_pts,2);
    
    head_time=find(head_pts(:,1) & ~same_pt);
    head_pts_sub=head_pts(head_time,:);
    framelist=1:max(head_time);
    head_pt_list=interp1(head_time,head_pts_sub,framelist,'pchip');
    
    tail_time=find(tail_pts(:,1) & ~same_pt);
    tail_pts_sub=tail_pts(tail_time,:); 
    tail_pt_list=interp1(tail_time,tail_pts_sub,framelist,'pchip');
  tips.head_pts=head_pt_list;
  tips.tail_pts=tail_pt_list;

 function C=processImage(handles)
 %for processing background subtraction and testing out different
     %filtering/processing methods     
     frameNumber=get(handles.slider1,'Value');
frameNumber=max(1,round(frameNumber));
Fid=getappdata(handles.figure1,'Fid');
    C=read(Fid,frameNumber);
    C=C(:,:,1);


     
    
background=getappdata(handles.figure1,'background');
frame_bg_lvl=getappdata(handles.figure1,'frame_bg_lvl');

background_raw=background(:,:,frame_bg_lvl(frameNumber));
    bg=normalizeRange(background_raw);
    bg_mask=bg>.5;
    bg_mask=AreaFilter(bg_mask,5000,[],8);
    bg_mask=imclose(bg_mask,true(12));
    bg_mask=imdilate(bg_mask,true(25));
    cc=bwconncomp(bg_mask);
    frame_bg_lvl(frameNumber)
    C=double(C);
    c=sum(sum(C.*background_raw))/sum(background_raw(:).^2);

  k=fspecial('Gaussian',5,2.5);
C2=imfilter(C,k);
[H,D]=hessianMatrix(C2,3);
[Heig,HeigVec]=hessianEig(H);

    C=double(C)-c*background_raw;
    
    se=strel('disk',11);
se=se.Neighborhood;

H1=stdfilt(HeigVec{1,1}.*Heig(:,:,1),se);
H2=stdfilt(HeigVec{2,1}.*Heig(:,:,1),se);
H=sqrt(H1.^2+H2.^2);

H=H*1000;
H=H-15;
H(H<0)=0;
C(C<0)=0;

%C=.3*H+.7*C;
C(C<0)=0;
C=smooth2a(C,5,5);


% --- Executes on button press in save_CL.
function save_CL_Callback(hObject, eventdata, handles)
% hObject    handle to save_CL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

eigenWormFile='eigenWorms.mat';
load(eigenWormFile);


centerline = getappdata(handles.figure1, 'centerline');
currentFolder=get(handles.currentFolder,'String');
dataFolder=fileparts(fileparts(currentFolder));
centerline=centerline{1};

wormcentered=FindWormCentered(centerline);
%project onto eigen basis

eigbasis=imresize(eigbasis,[size(eigbasis,1),size(wormcentered,1)]);


eigenProj=eigbasis*wormcentered;
%save outputs into behaviorAnalysis folder
behaviorFolder=[dataFolder filesep 'BehaviorAnalysis'];
display(['Saving data in ' behaviorFolder])
save([behaviorFolder filesep 'centerline'] ,'centerline','eigenProj'...
    ,'wormcentered');

save([behaviorFolder filesep 'CL_copy'] ,'centerline','eigenProj'...
    ,'wormcentered');


% --- Executes on button press in copyCL.
function copyCL_Callback(hObject, eventdata, handles)
% hObject    handle to copyCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i=1:str2double(get(handles.stepSize,'String'))
    frameNumber=get(handles.slider1,'Value');
    frameNumber=max(1,round(frameNumber));
    % use last frames CL
    CL_set=getappdata(handles.figure1,'centerline');
    cl = CL_set{1}(:,[2,1],frameNumber-1);
    %set current frame to last frames data
    CL_set{1}(:,[2,1],frameNumber)=cl;
    setappdata(handles.figure1,'centerline',CL_set);
    hold(handles.axes1,'on')
    plot(handles.axes1,cl(:,1),cl(:,2))
    hold(handles.axes1,'off');
    showImage(handles.slider1,eventdata)
    handles.slider1.Value=handles.slider1.Value+1;
end
     

