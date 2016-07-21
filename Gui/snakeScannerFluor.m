function varargout = snakeScannerFluor(varargin)
% SNAKESCANNERFLUOR MATLAB code for snakeScannerFluor.fig
%      SNAKESCANNERFLUOR, by itself, creates a new SNAKESCANNERFLUOR or raises the existing
%      singleton*.
%
%      H = SNAKESCANNERFLUOR returns the handle to a new SNAKESCANNERFLUOR or the handle to
%      the existing singleton*.
%
%      SNAKESCANNERFLUOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNAKESCANNERFLUOR.M with the given input arguments.
%
%      SNAKESCANNERFLUOR('Property','Value',...) creates a new SNAKESCANNERFLUOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before snakeScannerFluor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to snakeScannerFluor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help snakeScannerFluor

% Last Modified by GUIDE v2.5 22-Sep-2014 15:18:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @snakeScannerFluor_OpeningFcn, ...
                   'gui_OutputFcn',  @snakeScannerFluor_OutputFcn, ...
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


% --- Executes just before snakeScannerFluor is made visible.
function snakeScannerFluor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to snakeScannerFluor (see VARARGIN)

% Choose default command line output for snakeScannerFluor
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @showImage);


setappdata(handles.figure1,'holdaxes',false);
setappdata(handles.slider1,'hlistener',hlistener);
set(handles.slider1,'Max',2000);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes snakeScannerFluor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = snakeScannerFluor_OutputFcn(hObject, eventdata, handles) 
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

try
currentData=uipickfiles('FilterSpec',mostRecent);
catch
    currentData=uipickfiles();
end
currentData=currentData{1};

setappdata(0,'mostRecent',fileparts(currentData));
set(handles.currentFolder,'String',currentData);
aviFiles=dir([currentData filesep '*.avi']);
aviFiles={aviFiles.name}';
aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
if length(aviFiles)==2
aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
behaviorMovie=[currentData filesep aviFiles{~aviFluorIdx}];
fluorMovie=[currentData filesep aviFiles{aviFluorIdx}];
else
    display('Select avi files, behavior and then low mag fluor');
    movies=uipickfiles;
    behaviorMovie=movies{1};
    fluorMovie=movies{2};
end



[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(currentData);
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);


fluor2bfTimeOffset=fluorFlashTime-bfFlashTime;
fluorAll.frameTime=fluorAll.frameTime-fluor2bfTimeOffset(1);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

% make lookup tables for indices
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
%bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,bfAll.frameTime,'linear');
% 
% [hiImageIdx,ib]=unique(hiResData.imageIdx);
% hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));
% 
firstFullFrame=1;
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));


behaviorVidObj = VideoReader(behaviorMovie);
 fluorVidObj= VideoReader(fluorMovie);

 nFrames= behaviorVidObj.NumberOfFrames;
 
 display('Select Low Res Alignment')
lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
 
setappdata(handles.figure1,'nFrames',nFrames);
set(handles.slider1,'Value',firstFullFrame);
set(handles.slider1,'Max',nFrames);
set(handles.slider1,'Min',firstFullFrame);
set(handles.maxSlider,'String',num2str(nFrames));
set(handles.minSlider,'String','1');
setappdata(handles.figure1,'MaxValue',nFrames);
setappdata(handles.figure1,'vidObj',behaviorVidObj);
setappdata(handles.figure1,'fluorvidObj',fluorVidObj);
setappdata(handles.figure1,'lowResFluor2BF',lowResFluor2BF);
setappdata(handles.figure1,'fluorIdxLookup',fluorIdxLookup);
showImage(hObject)

% --- Executes on button press in flashSelect.
function flashSelect_Callback(hObject, eventdata, handles)
% hObject    handle to flashSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');

flashFile=uipickfiles('FilterSpec',mostRecent);

flashFile=flashFile{1};
imFlash=load(flashFile);
imFlash=imFlash.imFlash;
imFlash=imFlash-min(imFlash);
imFlash=imFlash>(max(imFlash)/2);
setappdata(handles.figure1,'imFlash',imFlash);



% --- Executes on button press in centerlineSelect.
function centerlineSelect_Callback(hObject, eventdata, handles)
% hObject    handle to centerlineSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');

centerlineFile=uipickfiles('FilterSpec',mostRecent);

centerlineFile=centerlineFile{1};
setappdata(handles.figure1,'centerlineFile',centerlineFile);
centerline=load(centerlineFile);
centerline=centerline.centerline;
setappdata(handles.figure1,'centerline',centerline);
set(handles.currentCenterline,'String',centerlineFile);

% --- Executes on button press in selectBackground.
function selectBackground_Callback(hObject, eventdata, handles)
% hObject    handle to selectBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mostRecent=getappdata(0,'mostRecent');

backgroundFile=uipickfiles('FilterSpec',mostRecent);
backgroundIm=imread(backgroundFile{1},'tiff');
backgroundIm=backgroundIm(:,:,1);
backgroundIm=double(backgroundIm)/256;

backgroundBW=backgroundIm>graythresh(backgroundIm);
bgIlevel=median(backgroundIm(backgroundBW));

setappdata(handles.figure1,'backgroundIm',backgroundIm);
setappdata(handles.figure1,'bacgroundBW',backgroundBW);
setappdata(handles.figure1,'bgIlevel',bgIlevel);




function showImage(hObject,eventdata)
handles=guidata(get(hObject,'Parent'));

backgroundIm=getappdata(handles.figure1,'backgroundIm');
bgIlevel=getappdata(handles.figure1,'bgIlevel');
backgroundBW=getappdata(handles.figure1,'backgroundBW');
winSize=50;
vidObj=getappdata(handles.figure1,'vidObj');
centerline=getappdata(handles.figure1,'centerline');
frameNumber=get(handles.slider1,'Value');
frameNumber=ceil(frameNumber);

   fluorVidObj=getappdata(handles.figure1,'fluorvidObj');
lowResFluor2BF=getappdata(handles.figure1,'lowResFluor2BF');
fluorIdxLookup=getappdata(handles.figure1,'fluorIdxLookup');
fluorFrameNumber=fluorIdxLookup(frameNumber);
fluorFrame = im2double((read(fluorVidObj, fluorFrameNumber)));
fluorFrame = fluorFrame(:,:,1);
  fluorFrame2=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
    'OutputView',lowResFluor2BF.Rsegment);



stretchSize=40;
I1 = im2double((read(vidObj, frameNumber)));
I1 = I1(:,:,1);
    bgscale=median(I1(backgroundBW))/bgIlevel;
bgscale=min(.9,bgscale);
I1=I1-backgroundIm*bgscale;
I1(I1<0)=0;
I1=imtophat(I1,strel('disk',30));
switch get(handles.channel,'Value')
    case 1
h=imagesc(I1,'Parent',handles.axes1);
    case 2
        h=imagesc(fluorFrame2,'Parent',handles.axes1);
end

%colormap(handles.axes1,'jet');
set(handles.currentFrame,'String',num2str(frameNumber));
CLcurrent=centerline(:,:,frameNumber);

CLcurrent=[interp1(CLcurrent(:,1),-stretchSize+1:100+stretchSize,'*linear','extrap')',...
    interp1(CLcurrent(:,2),-stretchSize+1:100+stretchSize,'*linear','extrap')'];

try
    hold(handles.axes1,'on');
    plot(handles.axes1,CLcurrent(:,2),CLcurrent(:,1),'black');
    hold(handles.axes1,'off');
catch
end

if get(handles.wormCoord,'Value')
          [lowResFluorInterp,newY,newX]=wormStraightening(CLcurrent,fluorFrame2,winSize,.25);
lowResFluorInterp=lowResFluorInterp';
       % lowResFluorInterp=interp2(fluorFrame2,newY,newX)';
       lowResFluorInterp=pedistalSubtract(lowResFluorInterp);
       lowResFluorInterp(isnan(lowResFluorInterp))=0;
    %   lowResFluorInterp=lowResFluorInterp.*(lowResFluorInterp>graythresh(lowResFluorInterp)/2);

    [Iinterp,newY,newX]=wormStraightening(CLcurrent,I1,winSize,.25);
    Iinterp=smooth2a(Iinterp',8,0);
    Iinterp=pedistalSubtract(Iinterp);
    Iinterp=Iinterp+9*lowResFluorInterp;
    lineWeight=nanmean(Iinterp);
    [~,Xgrid]=meshgrid(1:size(Iinterp,2),1:size(Iinterp,1));
    maxPos=bsxfun(@rdivide,sum(Xgrid.*Iinterp),sum(Iinterp));
%     [maxI,maxPos]=nanmax(Iinterp);
%     maxPos=((maxPos.*maxI)+winSize*nanmean(lineWeight)/2)./(maxI+nanmean(lineWeight)/2);
%maxPos=smooth(maxPos,5)';


newPos=[interp2(newY',1:length(maxPos),maxPos)' interp2(newX',1:length(maxPos),maxPos)'];

lowResFluorInterpThresh=(lowResFluorInterp>graythresh(lowResFluorInterp));
Isum=smooth(sum(lowResFluorInterp.*lowResFluorInterpThresh),3);
[ymax,imax,ymin,imin]=extrema(Isum);
%lowResFluorInterpThresh(:,imin)=0;
fluorCM=sum((Xgrid.*lowResFluorInterpThresh))./sum(lowResFluorInterpThresh);
[imax,ib]=sort(imax);
ymax=ymax(ib);
centroids=[imax fluorCM(imax)' ];

CL2=CLcurrent;

% stats=regionprops(lowResFluorInterpThresh,lowResFluorInterp,'Area','Centroid'...
%     ,'WeightedCentroid','WeightedCentroid');
% centroids=cell2mat({stats.Centroid}');
%centroids=sort(centroids,1);
%centroids=centroids([1,end],:);
maxPos=imax(ymax==max(ymax));

if maxPos>length(lowResFluorInterp)/2;
        centroids(imax<maxPos,:)=[];
%    centroids(centroids(:,1)<(maxPos-12),:)=[];
    centroids2=[interp2(newY',centroids(:,1),centroids(:,2)) ...
    interp2(newX',centroids(:,1),centroids(:,2))];
     CL2= CL2(1:round((maxPos-15)/4),:);
      CL2=[CL2;centroids2];

else
    centroids(imax>maxPos,:)=[];
     %   centroids(centroids(:,1)>(maxPos+12),:)=[];
        centroids2=[interp2(newY',centroids(:,1),centroids(:,2)) ...
    interp2(newX',centroids(:,1),centroids(:,2))];
  CL2=CL2(round((maxPos+15)/4):end,:);
  CL2=[centroids2;CL2];
end

CLdistance=[0;cumsum(sqrt(sum(diff(CL2).^2,2)))];
totLengthPix=max(CLdistance);
CL2=[interp1(CLdistance,CL2(:,1),1:1:totLengthPix)',...
    interp1(CLdistance,CL2(:,2),1:1:totLengthPix)'];
stretchSize=10;
CL2=[interp1(CL2(:,1),-stretchSize+1:length(CL2)+stretchSize,'*linear','extrap')',...
    interp1(CL2(:,2),-stretchSize+1:length(CL2)+stretchSize,'*linear','extrap')'];

[fluorFrameInterp2,newY,newX]=wormStraightening(CL2,fluorFrame2,winSize,.25);
CL2=[interp1(CL2(:,1),linspace(1,length(CL2),100),'*linear','extrap')',...
    interp1(CL2(:,2),linspace(1,length(CL2),100),'*linear','extrap')'];
    
length(CL2)
h2=imagesc(fluorFrameInterp2','Parent',handles.axes2);

       
       imagesc(lowResFluorInterpThresh.*lowResFluorInterp,'parent',handles.axes3);
             hold(handles.axes3,'on');
scatter(handles.axes3,centroids(:,1),centroids(:,2),'red')
             hold(handles.axes3,'off');
    
    hold(handles.axes2,'off')
    
   hold(handles.axes1,'on');
    plot(handles.axes1,CL2(:,2),CL2(:,1),'green');
    scatter(handles.axes1,centroids2(:,2),centroids2(:,1),'red')
    hold(handles.axes1,'off');
    
    

%   [CLnewY,CLnewX]=transformPointsInverse(lowResFluor2BF.t_concord...
%        , CLcurrent(:,2),CLcurrent(:,1));
%   [newYfluor,newXFluor]=transformPointsInverse(lowResFluor2BF.t_concord,...
%       newY,newX);

    
end






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
set(handles.slider1,'Value',currentFrame-1);
showImage(hObject)


% --- Executes on button press in forward1.
function forward1_Callback(hObject, eventdata, handles)
% hObject    handle to forward1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFrame=get(handles.slider1,'Value');
set(handles.slider1,'Value',currentFrame+1);
showImage(hObject)




% --- Executes on button press in snakeIt.
function snakeIt_Callback(hObject, eventdata, handles)
% hObject    handle to snakeIt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentFrame=ceil(get(handles.slider1,'Value'));
backgroundIm=getappdata(handles.figure1,'backgroundIm');
bgIlevel=getappdata(handles.figure1,'bgIlevel');
backgroundBW=getappdata(handles.figure1,'backgroundBW');
centerline=getappdata(handles.figure1,'centerline');
imFlash=getappdata(handles.figure1,'imFlash');
vidObj=getappdata(handles.figure1,'vidObj');
if get(handles.snakeIt,'Value')
    set(handles.snakeIt,'String','SNAKING')
else
    set(handles.snakeIt,'String','Snake It');
end

while get(handles.snakeIt,'Value')
    j=currentFrame;
if ~imFlash(j)
    I =im2double((read(vidObj, j)));
    I1 = I(:,:,1);
    bgscale=median(I1(backgroundBW))/bgIlevel;
bgscale=min(.9,bgscale);
I1=I1-backgroundIm*bgscale;
I1(I1<0)=0;

    I5 = im2double(im2double((read(vidObj, j+5))));
    I5 = I5(:,:,1);
    bgscale=median(I5(backgroundBW))/bgIlevel;
bgscale=min(.9,bgscale);
I5=I5-backgroundIm*bgscale;
I5(I5<0)=0;
% I1(I1(:)>quantile(I1(:),.999))=quantile(I1(:),.99);
% I5(I5(:)>quantile(I5(:),.99))=quantile(I5(:),.99);
I1=I1/std(I1(:));
I1(I1>10)=10;
I5=I5/std(I5(:));
I5(I5>10)=10;
    I1=normalizeRange(I1);
    I5=normalizeRange(I5);
    [centerline(:,:,j),  Icirc, realnumtips(1,j)] = bettersnakeit_jn(centerline(:,:,j-1), I1, I5, get(handles.manualTip,'Value'));
    circAll(:,:,j)=Icirc;
    CL = centerline(:,:,j);
  

    Icomb = (Icirc)+(I1);
else
    centerline(:,:,j)=centerline(:,:,j-1);
    realnumtips(1,j)=realnumtips(1,j-1);
end

    
    imagesc(Icomb, 'parent',handles.axes1);
    hold(handles.axes1,'on');
    plot(handles.axes1,centerline(:,2,j), centerline(:,1, j), '-g')
    plot(handles.axes1,centerline(1,2,j), centerline(1,1,j), 'og')
    plot(handles.axes1,centerline(100,2,j), centerline(100,1,j), 'ob')
    hold(handles.axes1,'off');
   % text(10, 10, num2str((j)), 'Color', 'y');
        %lnot = length;
        currentFrame=currentFrame+1;
        set(handles.slider1,'Value',currentFrame);
end

if get(handles.snakeIt,'Value')
    set(handles.snakeIt,'String','SNAKING')
else
    set(handles.snakeIt,'String','Snake It');
end
setappdata(handles.figure1,'centerline',centerline);


% --- Executes on key press with focus on snakeIt and none of its controls.
function snakeIt_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to snakeIt (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in manualTip.
function manualTip_Callback(hObject, eventdata, handles)
% hObject    handle to manualTip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manualTip


% --- Executes on button press in saveCenterline.
function saveCenterline_Callback(hObject, eventdata, handles)
% hObject    handle to saveCenterline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

centerlineFile=get(handles.currentCenterline,'String');
centerline=getappdata(handles.figure1,'centerline');
try
save(centerlineFile,'centerline','-append');
catch
    save(centerlineFile,'centerline');

    
end

    


function currentCenterline_Callback(hObject, eventdata, handles)
% hObject    handle to currentCenterline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentCenterline as text
%        str2double(get(hObject,'String')) returns contents of currentCenterline as a double


% --- Executes during object creation, after setting all properties.
function currentCenterline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentCenterline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in wormCoord.
function wormCoord_Callback(hObject, eventdata, handles)
% hObject    handle to wormCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wormCoord


% --- Executes on selection change in channel.
function channel_Callback(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel


% --- Executes during object creation, after setting all properties.
function channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
