function [dat_images,avi0_images,avi1_images]=ExtractAlignmentImagesFromVideo(imFolder)

if nargin==0
    mostRecent=getappdata(0,'mostRecent');
    imFolder=uipickfiles('filterspec',mostRecent);
    imFolder=imFolder{1};
    setappdata(0,'mostRecent',imFolder);
end
%% cut up dat file if present
h=dir([imFolder filesep '*.dat']);
if ~isempty(h)
    datFile=[imFolder filesep h(1).name];
    [rows, cols]=getdatdimensions(datFile);
    
    camData=importdata([imFolder filesep 'CameraFrameData.txt']);
    labJackData=importdata([imFolder filesep 'LabJackData.txt']);
    frameno=camData.data(:,1);
    savedFrames=camData.data(:,2);
    time=frameno(diff(savedFrames)>0);
    
    timestep=median(diff(time));
    frame_breaks=diff(time)>(timestep*10);
    frame_index=cumsum(frame_breaks)+1;
    
    Fid=fopen(datFile);
    status=fseek(Fid,0,1);
    stackSize=floor(ftell(Fid)/(2*rows*cols)-1);
    status=fseek(Fid,0,-1);
    
    dat_images=zeros(rows,cols,max(frame_index));
    progressbar(0);
    for idx=1:max(frame_index)
    progressbar(idx/max(frame_index));
        display(['starting index ' num2str(idx)])
        
        fileNamedat=[imFolder filesep 'dat_' num2str(idx) '.tif'];
        if exist(fileNamedat,'file')
            medImageDat=importdata(fileNamedat);
        else
            range=[find(frame_index==idx,1,'first'),...
                find(frame_index==idx,1,'last')];
            rowSearch=diff(range);
            pixelValues=fread(Fid,rows*cols*rowSearch,'uint16',0,'l');
            pixelValues=reshape(pixelValues,rows,cols,rowSearch);
            medImageDat=median(pixelValues,3);
            tiffwrite(fileNamedat,single(medImageDat),'tif');
        end
        dat_images(:,:,idx)=medImageDat;
    end
    
    fclose(Fid);
end

%% do the same if avi is present
h=dir([imFolder filesep '*.avi']);
if isempty(h)
    lowMagFolder=dir([imFolder filesep 'LowMag*']);
    if ~isempty(lowMagFolder)
        imFolder=[imFolder filesep lowMagFolder(1).name];
    else
        error('No AVI files found!')
    end
end

camData=importdata([imFolder filesep 'CamData.txt']);
time=camData.data(:,2);
timestep=median(diff(time));
frame_breaks=diff(time)>(timestep*10);
frame_index=cumsum(frame_breaks)+1;

cam0File=[imFolder filesep 'cam0.avi'];
cam1File=[imFolder filesep 'cam1.avi'];

vidObj0 = VideoReader(cam0File);
vidObj1 = VideoReader(cam1File);

avi0_images=zeros(vidObj0.Width,vidObj0.Height,max(frame_index));
avi1_images=zeros(vidObj1.Width,vidObj1.Height,max(frame_index));
    progressbar(0);

for idx=1:max(frame_index)
    progressbar(idx/max(frame_index))
    display(['starting index ' num2str(idx)])

    fileName0=[imFolder filesep 'cam0_' num2str(idx)];
    fileName1=[imFolder filesep 'cam1_' num2str(idx)];
    
    %if files exist, read them, if they dont, write them
    if exist(fileName0,'file') && exist(fileName1,'file')
        meanImage0=imread(fileName0);
        meanImage1=imread(fileName1);
    else
        range=[find(frame_index==idx,1,'first'),...
            find(frame_index==idx,1,'last')];
        temp0= read(vidObj0, range);
        temp1= read(vidObj1, range);
        
        meanImage0=squeeze(mean(temp0,4));
        meanImage1=squeeze(mean(temp1,4));
        tiffwrite(fileName0,single(meanImage0),'tif');
        tiffwrite(fileName1,single(meanImage1),'tif');
    end
    avi0_images(:,:,idx)=meanImage0;
    avi1_images(:,:,idx)=meanImage1;
    
end


