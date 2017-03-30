imFolder=uipickfiles();
imFolder=imFolder{1};

%% cut up dat file if present
h=dir([imFolder filesep '*.dat']);
if ~isempty(h)
    datFile=[imFolder filesep h(1).name];
[rows, cols]=getdatdimensions(datFile);
    
camData=importdata([imFolder filesep 'cameraFrameData.txt']);
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

for idx=1:max(frame_index)
    display(['starting index ' num2str(idx)])
    range=[find(frame_index==idx,1,'first'),...
        find(frame_index==idx,1,'last')];
    rowSearch=diff(range);
pixelValues=fread(Fid,rows*cols*rowSearch,'uint16',0,'l');
  pixelValues=reshape(pixelValues,rows,cols,rowSearch);
  medImageDat=median(pixelValues,3);
        fileNamedat=[imFolder filesep 'dat_' num2str(idx) '.tif'];
tiffwrite(fileNamedat,single(medImageDat),'tif');

end





fclose(Fid);
end

    %% do the same if avi is present
    h=dir([imFolder filesep '*.avi']);
if ~isempty(h)
camData=importdata([imFolder filesep 'camData.txt']);
time=camData.data(:,2);
timestep=median(diff(time));
frame_breaks=diff(time)>(timestep*10);
frame_index=cumsum(frame_breaks)+1;

cam0File=[imFolder filesep 'cam0.avi'];
cam1File=[imFolder filesep 'cam1.avi'];

        vidObj0 = VideoReader(cam0File);
        vidObj1 = VideoReader(cam1File);

for idx=1:max(frame_index)
    display(['starting index ' num2str(idx)])
    range=[find(frame_index==idx,1,'first'),...
        find(frame_index==idx,1,'last')];
        temp0= read(vidObj0, range);
        temp1= read(vidObj1, range);
        
        meanImage0=squeeze(mean(temp0,4));
        meanImage1=squeeze(mean(temp1,4));
        fileName0=[imFolder filesep 'cam0_' num2str(idx)];
        fileName1=[imFolder filesep 'cam1_' num2str(idx)];
tiffwrite(fileName0,single(meanImage0),'tif');
tiffwrite(fileName1,single(meanImage1),'tif');

end

end
