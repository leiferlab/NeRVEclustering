% spt_interactive is used to determine the proper settings for located
% fluorescent features in a single particle image or movie.  It will
% display the image(s) with positions at both the pixel and sub-pixel level
% accuracies.  Parameters may be changed in order to view their affect on
% the detected positions
function varargout = spt_interactive(varargin)
% INPUTS:
% sptpara = optional input containing spt parameters *** Need to add
%           exact values required***

% OUTPUTS:
% sptpara = final spt parameters explored *** Need to add exact values in
%           output
%
% Dependencies: tiffread.m
%               bpass.m
%               pkfnd.m
%               centfind.m
%
% Author: Colin Ingram
% Created: Nov 2009
% Version: 1.00
%
% Revisions:
% Version     Date  Author   Description


% spt_interactive M-file for spt_interactive.fig
%      spt_interactive, by itself, creates a new spt_interactive or raises the existing
%      singleton*.
%
%      H = spt_interactive returns the handle to a new spt_interactive or the handle to
%      the existing singleton*.
%
%      spt_interactive('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in spt_interactive.M with the given input arguments.
%
%      spt_interactive('Property','Value',...) creates a new spt_interactive or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spt_interactive_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spt_interactive_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spt_interactive

% Last Modified by GUIDE v2.5 15-Apr-2010 14:26:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spt_interactive_OpeningFcn, ...
                   'gui_OutputFcn',  @spt_interactive_OutputFcn, ...
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



% Execution and Annihilation
% Executes just before spt_interactive is made visible.
function spt_interactive_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spt_interactive (see VARARGIN)

% Choose default command line output for spt_interactive
handles.output = hObject;

% create play mode timer
playt.TimerFcn = {@TmrFcn,handles.figure1};
playt.BusyMode = 'Queue';
playt.ExecutionMode = 'FixedRate';
playt.Period = 1/str2double(get(handles.frame_rate,'String'));
% playt.StopFcn = {@TmrStop,handles.figure1};
handles.tmr = timer(playt);

% Handle inputs
if nargin == 3
    
    % spt_interactive called without inputs request filenames
    [filename,filepath]=uigetfile({'*.tif', 'Tiff File (*.tif)';...
        '*.stk', 'Metamorph Stack (*.stk)'},...
        'Choose Movie(s) to Analyze',...
        'MultiSelect','on');
    cd(filepath);
    
    if iscell(filename)
        ntrials = size(filename,2);
        stackn = cell(1,ntrials);
        for i = 1:ntrials
            stackn{i} = filename{i};
        end
        
    else
%         ntrials = 1;
        stackn{1} = filename;
    end
    
    % get settings
    prompt = {'Particle Threshold:',...
        'Integrated Threshold[optional]:',...
        'Particle Diameter(pixel) [must be odd]:',...
        'Center Finding Method [0 for 2D Gaussian, 1 for centroid]',...
        ['Windows Size (pixel) [Radius of sub-pixel determination',...
        'window size, Particle Diameter + 4 is a good value.',...
        'Must be odd!']};
    u_name = 'Input Parameters for Interactive Feature Location';
    numlines = 1;
    defaultanswer = {'','0','3','1','7'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    user_var = inputdlg(prompt,u_name,numlines,defaultanswer,options);
     
    % setup sptpara
    sptpara.thresh = str2double(user_var{1});
    sptpara.IntTh = str2double(user_var{2});
    sptpara.dia = str2double(user_var{3});
    sptpara.fitmethod = str2double(user_var{4});
    sptpara.boxr = str2double(user_var{5});
    sptpara.intfile = stackn;
    
    % error handling on settings
    while true
        if ~rem(sptpara.dia,2)
            
            % Warning Box
            warnstr = 'Particle Diameter must be odd! Please Reset!';
            dlgname = 'SPT Warning!';
            createmode.Interpreter = 'tex';
            createmode.WindowStyle = 'model';
            h = warndlg(warnstr,dlgname,createmode);
            uiwait(h);
            
            % new input
            prompt = 'Particle Diameter(pixel, must be odd):';
            u_name = 'Input Parameters for Interactive Feature Location';
            numlines = 1;
            defaultanswer = {num2str(sptpara.dia)};
            options.Resize = 'on';
            options.WindowStyle = 'normal';
            options.Interpreter = 'tex';
            user_var = inputdlg(prompt,u_name,numlines,defaultanswer,options);
            
            sptpara.dia = str2double(user_var{1});
        else
            break
        end
    end
    
    while true
        if ~rem(sptpara.boxr,2)
            
            % Warning Box
            warnstr = 'Windows Size must be odd! Please Reset!';
            dlgname = 'SPT Warning!';
            createmode.Interpreter = 'tex';
            createmode.WindowStyle = 'model';
            h = warndlg(warnstr,dlgname,createmode);
            uiwait(h);
            
            % new input
            prompt = 'Windows Size(pixel, must be odd):';
            u_name = 'Input Parameters for Interactive Feature Location';
            numlines = 1;
            defaultanswer = {num2str(sptpara.boxr)};
            options.Resize = 'on';
            options.WindowStyle = 'normal';
            options.Interpreter = 'tex';
            user_var = inputdlg(prompt,u_name,numlines,defaultanswer,options);
            
            sptpara.boxr = str2double(user_var{1});
        else
            break
        end
    end

    
elseif nargin > 3    
    
    sptpara = varargin{1};

end

thresh = sptpara.thresh;
IntTh = sptpara.IntTh;
dia = sptpara.dia;
opt.method = sptpara.fitmethod;
boxr = sptpara.boxr;

% Settings Panel
set(handles.files,'String',sptpara.intfile);
set(handles.threshold,'String',num2str(thresh),'UserData',num2str(thresh));
set(handles.IntTh,'String',num2str(IntTh),'UserData',num2str(IntTh));
set(handles.diameter,'String',num2str(dia),'UserData',num2str(dia));
set(handles.box,'String',num2str(boxr),'UserData',num2str(boxr));
set(handles.method,'String',num2str(opt.method),'UserData',num2str(opt.method));
set(handles.setpan,'UserData',sptpara)
handles.SeteditH = [handles.threshold,handles.IntTh,handles.diameter,handles.box,handles.method];
handles.SetbuttH = [handles.load,handles.resetPrev,handles.resetOrg];

% load image
handles.filename = sptpara.intfile{1};

% Get information about the tiff file
info = imfinfo(handles.filename);
imcnt = numel(info);
handles.imcnt = imcnt;

% read the first plane to get general tags
handles.im = tiffread(handles.filename,1,1);

% preallocate space for image
handles.im.rawImg = zeros([info(1).Height,info(1).Width,imcnt],...
    class(handles.im.data));

% remove data field
handles.im = rmfield(handles.im,'data');

% read in each plane
for jIFD = 1:imcnt
    handles.im.rawImg(:,:,jIFD) = imread(handles.filename, jIFD,...
        'Info', info);
end

% initial data
[handles.im,handles.r] = calc_data(handles);

% Plot
setappdata(handles.axes1,'first_call',0)
handles.frame=1;
handles = plot_frame(handles);

% update data
guidata(hObject, handles);


% Outputs from this function are returned to the command line.
function varargout = spt_interactive_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% wait for close to excute Output Function
uiwait(handles.figure1);

% output set in exit function
handles = guidata(hObject);
varargout{1} = handles.output;
delete(handles.tmr);
delete(hObject);

% Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: delete(hObject) closes the figure
sptpara = get(handles.setpan,'UserData');
sptpara.thresh = str2double(get(handles.threshold,'String'));
sptpara.IntTh = str2double(get(handles.IntTh,'String'));
sptpara.dia = str2double(get(handles.diameter,'String'));
sptpara.boxr = str2double(get(handles.box,'String'));
sptpara.fitmethod = str2double(get(handles.method,'String'));
handles.output = sptpara;
guidata(hObject, handles);
uiresume
%delete(hObject);

% Executes on button press in ExitButton.
function exitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(handles.figure1,eventdata,handles);



% Plot Tools
% Executes on PosToggle.
function PosToggle_Callback(hObject, eventdata, handles)
% hObject    handle to PosToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % Toggle button is pressed, take appropriate action
    set(hObject,'String','Show Positions');
    set(hObject,'ForegroundColor',[0 0 1]);
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed, take appropriate action
    set(hObject,'String','Hide Positions');
    set(hObject,'ForegroundColor',[1 0 0]);
end

% Update handles structure
guidata(hObject, handles);

% Plot
plot_frame(handles);

% Executes on RawToggle.
function RawToggle_Callback(hObject, eventdata, handles)
% hObject    handle to RawToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');

% turn on raw
if button_state == get(hObject,'Max')
    % Toggle button is pressed, take appropriate action
    set(hObject,'String','Show Filtered Image');
    set(hObject,'ForegroundColor',[0 0 1]);    
    handles.Climfilter = get(handles.axes1,'Clim');
    set(handles.axes1,'Clim',handles.Climraw)
    
% turn off filtered
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed, take appropriate action
    set(hObject,'String','Show Raw Image');
    set(hObject,'ForegroundColor',[1 0 0]);
    handles.Climraw = get(handles.axes1,'Clim');
    set(handles.axes1,'Clim',handles.Climfilter)    
end

% Update handles structure
guidata(hObject, handles);

% Plot
plot_frame(handles);

% Executes on Advance Button.
function advance_Callback(hObject, eventdata, handles)
% hObject    handle to advance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.frame<handles.imcnt
    handles.frame = handles.frame+1;
    guidata(hObject, handles);
end

% Plot
plot_frame(handles);

% Executes on Previous Button.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.frame>1
    handles.frame = handles.frame-1;
    guidata(hObject, handles);
end

% Plot
plot_frame(handles);

% Helper Function for Plotting
function handles = plot_frame(handles)
% Inputs
pos_state = get(handles.PosToggle,'Value');
raw_state = get(handles.RawToggle,'Value');
Clim = get(handles.axes1,'Clim');

frame = handles.frame;
r = handles.r;
imcnt = handles.imcnt;
file = handles.filename;


% plot
% if raw button is pressed plot raw image
if raw_state == get(handles.RawToggle,'Min')
    imdata = handles.im.filterImg(:,:,frame);
elseif raw_state == get(handles.RawToggle,'Max')
    imdata = handles.im.rawImg(:,:,frame);
end

handles.image1 = imagesc(imdata);
colormap(gray);

% if plot_frame is being called for the first time for a particlar file
% set axis properties
if ~(getappdata(handles.axes1,'first_call') == 1)
    axis image
handles.Climraw = [min(min(handles.im.rawImg(:,:,1))),...
    max(max((handles.im.rawImg(:,:,1))))];
handles.Climfilter = [min(min(handles.im.filterImg(:,:,1))), ...
    max(max((handles.im.filterImg(:,:,1))))];
%     handles.Climraw = [min(handles.im(1).data(:)) max(handles.im(1).data(:))];
%     handles.Climfilter = [min(handles.im(1).filter(:)) max(handles.im(1).filter(:))];
elseif get(handles.autoscale,'Value')== get(handles.autoscale,'Min');
    set(handles.axes1,'Clim',Clim);
end

% if position button is not pressed plot positions
if pos_state == get(handles.PosToggle,'Min') && ~(isempty(handles.im.planeAttr(frame).particlepxl))

    % plot pixel threshold locations
    coursex = handles.im.planeAttr(frame).particlepxl(:,1);
    coursey = handles.im.planeAttr(frame).particlepxl(:,2);
%     pksz = size(handles.im.planeAttr(frame).particlepxl,1);
    hold on
    plot(coursex,coursey,'y.');
    
    % plot integrated threshold locations
    intcoursex = handles.im.planeAttr(frame).particleintthr(:,1);
    intcoursey = handles.im.planeAttr(frame).particleintthr(:,2);
    plot(intcoursex,intcoursey,'g.');
    
    
    % if boxr is much greater the dia sometime there can be no cents
    if ~(isempty(handles.im.planeAttr(frame).particlepos))
        finex = handles.im.planeAttr(frame).particlepos(:,1);
        finey = handles.im.planeAttr(frame).particlepos(:,2);
        plot(finex,finey,'rx');
        plot(finex,finey,'ro');
        
        % 4/15 changed box to finexy because with new threshold there can be
        % lots of coursexy points
        pksz = size(handles.im.planeAttr(frame).particlepos,1);
        for k = 1:pksz
            x=round(finex(k));
            y=round(finey(k));
            boxx = [x-r,x+r,x+r,x-r,x-r];
            boxy = [y+r,y+r,y-r,y-r,y+r];
            plot(boxx,boxy,'b-');
        end
    end
    hold off
end

% if cursor mode is on reinitialize
if strcmp(get(handles.cursortoggle,'State'),'on')
    cursortoggle_OnCallback(handles.cursortoggle, 1, handles);
    handles = guidata(handles.figure1);
end

set(handles.axes1,'NextPlot','replacechildren');
setappdata(handles.axes1,'first_call',1);
set(handles.ftxt,'String',['Frame ',num2str(frame),' of ',num2str(imcnt)]);
set(handles.txtFile,'String',file);

% Update handles structure
guidata(handles.figure1, handles);


% Settings Table and New Settings Updates
% threshold error handling
function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));
if value < 0
    warndlg({'Threshold must be larger than zero'},'Invalid Entry');
    set(hObject,'String',get(hObject,'UserData'));
end

% Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% IntTh Error Handling... there is none for now
function IntTh_Callback(hObject, eventdata, handles)
% hObject    handle to IntTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IntTh as text
%        str2double(get(hObject,'String')) returns contents of IntTh as a double


% --- Executes during object creation, after setting all properties.
function IntTh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% diameter box error handling
function diameter_Callback(hObject, eventdata, handles)
% hObject    handle to diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameter as text
%        str2double(get(hObject,'String')) returns contents of diameter as a double
value = str2double(get(hObject,'String'));
if ~rem(value,2)
    warndlg({'Particle Diameter must be Odd'},'Invalid Entry');
    set(hObject,'String',get(hObject,'UserData'));
end

% Executes during object creation, after setting all properties.
function diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% box size error handling
function box_Callback(hObject, eventdata, handles)
% hObject    handle to box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box as text
%        str2double(get(hObject,'String')) returns contents of box as a double
value = str2double(get(hObject,'String'));
if ~rem(value,2)
    warndlg({'Box Size must be Odd'},'Invalid Entry');
    set(hObject,'String',get(hObject,'UserData'));
end

% Executes during object creation, after setting all properties.
function box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Method Box error Handling
function method_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of method as text
%        str2double(get(hObject,'String')) returns contents of method as a double
value = str2double(get(hObject,'String'));
if ~(value == any([0,1]))
    warndlg({'Please enter 0 to use Gaussian or 1 to use Centroid'},'Invalid Entry');
    set(hObject,'String',get(hObject,'UserData'));
end

% Executes during object creation, after setting all properties.
function method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on Load New Settings to Caluclate new data with new settings
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = waitbar(0,'Calculated Positions With New Settings');
[handles.im,handles.r] = calc_data(handles);
close(h);

% Toggle edit button
set(handles.editSet,'String','Edit Settings');
set(handles.editSet,'BackgroundColor',[0 1 0]);
set(handles.editSet,'Value',get(handles.editSet,'Min'));

% disable txtbox
set(handles.SeteditH,'Enable','inactive');

% hide buttons
set(handles.SetbuttH,'Visible','off');

% enable plot tools
set(handles.PosToggle,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% Plot
plot_frame(handles);

% Executes on Reset Previous Button to reset previous Settings
function resetPrev_Callback(hObject, eventdata, handles)
% hObject    handle to resetPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.threshold,'String',get(handles.threshold,'UserData'));
set(handles.IntTh,'String',get(handles.IntTh,'UserData'));
set(handles.diameter,'String',get(handles.diameter,'UserData'));
set(handles.box,'String',get(handles.box,'UserData'));
set(handles.method,'String',get(handles.method,'UserData'));
    
% Executes resetOrg button to Reset Orginal Settings
function resetOrg_Callback(hObject, eventdata, handles)
% hObject    handle to resetOrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Handle inputs
sptpara = get(handles.setpan,'UserData');
thresh = sptpara.thresh;
IntTh = sptpara.IntTh;
dia = sptpara.dia;
opt.method = sptpara.fitmethod;
boxr = sptpara.boxr+sptpara.dia;

% Settings Panel
set(handles.threshold,'String',num2str(thresh));
set(handles.IntTh,'String',num2str(IntTh));
set(handles.diameter,'String',num2str(dia));
set(handles.box,'String',num2str(boxr));
set(handles.method,'String',num2str(opt.method));

% Executes editSet button to enable Settings Editting
function editSet_Callback(hObject, eventdata, handles)
% hObject    handle to editSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % Toggle button is pressed, take appropriate action
    set(hObject,'String','Cancel');
    set(hObject,'BackgroundColor',[1 0 0]);
    
    % enable txtbox
    set(handles.SeteditH,'Enable','on');
    set(handles.threshold,'UserData',get(handles.threshold,'String'));
    set(handles.IntTh,'UserData',get(handles.IntTh,'String'));
    set(handles.diameter,'UserData',get(handles.diameter,'String'));
    set(handles.box,'UserData',get(handles.box,'String'));
    set(handles.method,'UserData',get(handles.method,'String'));
    
    % show buttons
    set(handles.SetbuttH,'Visible','on');
    
    % disable pos button
    set(handles.PosToggle,'Enable','off');
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed, take appropriate action
    set(hObject,'String','Edit Settings');
    set(hObject,'BackgroundColor',[0 1 0]);
    
    % disable txtbox
    set(handles.SeteditH,'Enable','inactive');
    set(handles.threshold,'String',get(handles.threshold,'UserData'));
    set(handles.IntTh,'String',get(handles.IntTh,'UserData'));
    set(handles.diameter,'String',get(handles.diameter,'UserData'));
    set(handles.box,'String',get(handles.box,'UserData'));
    set(handles.method,'String',get(handles.method,'UserData'));
    
    % hide buttons
    set(handles.SetbuttH,'Visible','off');
    
    % enable plot tools
    set(handles.PosToggle,'Enable','on');
    
end

% Helper Function for Calculating Positions
function [im,r] = calc_data(handles)
sptpara.thresh = str2double(get(handles.threshold,'String'));
sptpara.IntTh = str2double(get(handles.IntTh,'String'));
sptpara.dia = str2double(get(handles.diameter,'String'));
sptpara.boxr = str2double(get(handles.box,'String'));
sptpara.fitmethod = str2double(get(handles.method,'String'));
sptpara.file = handles.filename;
r = sptpara.boxr/2;

% load image
im = handles.im;
% imcnt = handles.imcnt;


if isnan(sptpara.thresh+sptpara.IntTh+sptpara.dia+sptpara.boxr+sptpara.fitmethod)
    return
end
    
% calculate
[im, ~] = spt_fndpos(sptpara,im,0);
% bpb not sure why this section is commented out, it has not been corrected
% for new im structure
% for i=1:imcnt
%     im(i).filter = bpass(double(im(i).data),1,dia+2);
%     im(i).pkfind = pkfnd(im(i).filter,thresh,dia+2);
%     if isempty(im(i).pkfind)
%         im(i).cent = [];
%     else
%         im(i).cent = centfind(im(i).filter,im(i).pkfind,boxr,opt);
%     end
% end


% --- Executes on selection change in files.
function files_Callback(hObject, eventdata, handles)
% hObject    handle to files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files
ind = get(hObject,'Value');
list = get(hObject,'String');
handles.filename = list{ind};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in newfile.
function newfile_Callback(hObject, eventdata, handles)
% hObject    handle to newfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = handles.filename;

h = waitbar(0,'Loading New Image');
% load image

% Get information about the tiff file
info = imfinfo(handles.filename);
imcnt = numel(info);
handles.imcnt = imcnt;

% read the first plane to get general tags
handles.im = tiffread(handles.filename,1,1);

% preallocate space for image
handles.im.rawImg = zeros([info(1).Height,info(1).Width,imcnt],...
    class(handles.im.data));

% remove data field
handles.im = rmfield(handles.im,'data');

% read in each plane
for jIFD = 1:imcnt
    handles.im.rawImg(:,:,jIFD) = imread(handles.filename, jIFD,...
        'Info', info);
end



% initial data
[handles.im,handles.r] = calc_data(handles);
close(h);

% Plot
setappdata(handles.axes1,'first_call',0)
handles.frame=1;
plot_frame(handles);

% update data
guidata(hObject, handles);

%Timer Function
function TmrFcn(src,event,handles)
handles = guidata(handles);
frame = handles.frame;
imcnt = handles.imcnt;
loop = get(handles.loop,'Value');
if frame < imcnt
    handles.frame = frame + 1;
elseif loop == get(handles.loop,'Max')
    handles.frame = 1;
else
    set(handles.playToggle,'Value',get(handles.playToggle,'Min'));
    playToggle_Callback(handles.playToggle, event, handles);
end
guidata(handles.figure1,handles);
plot_frame(handles);

% %Timer Stop
% function TmrStop(src,event,handles) 
% handles = guidata(handles);
% set(handles.playToggle,'String','Play');
% set(hObject,playToggle,'ForegroundColor',[0 1 0]);

% --- Executes on button press in playToggle.
function playToggle_Callback(hObject, eventdata, handles)
% hObject    handle to playToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of playToggle
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % Toggle button is pressed, take appropriate action
    set(hObject,'String','Stop');
    set(hObject,'ForegroundColor',[1 0 0]);
    set(handles.cursortoggle,'State','off'); % having cursor on creates errs
    start(handles.tmr)
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed, take appropriate action
    set(hObject,'String','Play');
    set(hObject,'ForegroundColor',[0 1 0]);
    stop(handles.tmr)
    set(handles.cursortoggle,'State','on'); % cursor on is default!
end

% --- Executes on button press in writebutton.
function writebutton_Callback(hObject, eventdata, handles)
% hObject    handle to writebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loop.
function loop_Callback(hObject, eventdata, handles)
% hObject    handle to loop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loop



function frame_rate_Callback(hObject, eventdata, handles)
% hObject    handle to frame_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_rate as text
%        str2double(get(hObject,'String')) returns contents of frame_rate as a double
fps = str2double(get(hObject,'String'));
period = 1./fps;
set(handles.tmr,'Period',period);

% --- Executes during object creation, after setting all properties.
function frame_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in first.
function first_Callback(hObject, eventdata, handles)
% hObject    handle to first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.frame = 1;
plot_frame(handles);

% --- Executes on button press in last.
function last_Callback(hObject, eventdata, handles)
% hObject    handle to last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.frame = handles.imcnt;
plot_frame(handles);

% --- Executes on button press in GoTo.
function GoTo_Callback(hObject, eventdata, handles)
% hObject    handle to GoTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.frame = str2double(get(handles.GoToBox,'String'));
plot_frame(handles);

function GoToBox_Callback(hObject, eventdata, handles)
% hObject    handle to GoToBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GoToBox as text
%        str2double(get(hObject,'String')) returns contents of GoToBox as a
%        double

% protect them from themselves
frame = str2double(get(hObject,'String'));
if frame < 1 || isnan(frame)
    set(hObject,'String','1');
elseif frame > handles.imcnt
    set(hObject,'String',num2str(handles.imcnt));
else
    set(hObject,'String',round(frame));
end
    

% --- Executes during object creation, after setting all properties.
function GoToBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GoToBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function cursortoggle_OnCallback(hObject, eventdata, handles)
% hObject    handle to cursortoggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% SetUp Cursor Mode
handles.pxvalue = impixelinfoval(handles.figure1,handles.image1);
set(handles.pxvalue,'HorizontalAlignment','right')
set(handles.pxvalue,'FontSize',12)
set(handles.pxvalue,'FontWeight','bold')
set(handles.pxvalue,'ForegroundColor',[1 0 0])
set(handles.pxvalue,'Units','normalized')
% set(handles.pxvalue,'Position',[1200,3,139,23])
set(handles.pxvalue,'Position',[.8861,.0026,.1027,.0294])
guidata(hObject, handles)

% --------------------------------------------------------------------
function cursortoggle_OffCallback(hObject, eventdata, handles)
% hObject    handle to cursortoggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.pxvalue);


% --- Executes on button press in scalebutton.
function scalebutton_Callback(hObject, eventdata, handles)
% hObject    handle to scalebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=imcontrast;
uiwait(h);


% --- Executes on button press in autoscale.
function autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoscale
