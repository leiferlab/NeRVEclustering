function  imFlash=findDatFlash(datFile,rows,cols,rowSearch)
% function takes a .dat file for worm data and searches each of the images
% for a flash. The images must be rowxcol matrix of uint16, works super
% fast. Can also select input rows. This will find the flash signal by
% averaging over the first rows in order to speed up program. The flash is
% normally bright enough that I only need 1 row or less...

%% load images and find flash in images using user defined ROI
if nargin==0
    datFile=uipickfiles;
    datFile=datFile{1};
end

if nargin<3
    try
        [rows, cols]=getdatdimensions(datFile);
    catch
    rows=1200;
    cols=600;
    end
end
if nargin<4
    rowSearch=rows;
end




if strfind(datFile, '.dat');
    Fid=fopen(datFile);
status=fseek(Fid,0,1);
stackSize=floor(ftell(Fid)/(2*rows*cols));
status=fseek(Fid,0,-1);
else
    error('File must be a .dat binary file in uint16 form');
end

% if flag.customRoi
% progressbar;
% for i=round(stackSize/2):round(stackSize/2)+200;
%     if aviFlag;
%         temp= read(vidObj, i);
% temp=(sum(double(temp),3));
%     else
%     temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
%     end
%     initialIm=max(initialIm,temp);
%     progressbar((i-round(stackSize/2))/200);
% 
% end
% 
% imsize=size(initialIm);
% fig=imagesc(initialIm);
% display('Select area to find flash');
% roiFlash=roipoly;
% delete(fig)
% else
%     roiFlash=true(size(initialIm));
% end

% fig=imagesc(initialIm);
% display('Select an cell to try to align');
% roiCell1=roipoly;
% delete(fig)
% 
% fig=imagesc(initialIm);
% display('Select another');
% roiCell2=roipoly;
% close all


%% search for the flash
imFlash=zeros(1,stackSize);
%cell1=imFlash;
%cell2=imFlash;
progressbar(0)
for iFrame=1:stackSize
    progressbar(iFrame/stackSize)
    pixelValues=fread(Fid,rows*rowSearch,'uint16',0,'l');
    status=fseek(Fid,2*cols*(rows)*iFrame,-1);

imFlash(iFrame)=mean(pixelValues);
  %  status=fseek(Fid,2*frameNumber*1024^2,-1);

end

%%
fclose(Fid);
