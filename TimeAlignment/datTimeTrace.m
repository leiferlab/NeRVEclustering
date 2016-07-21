function  imFlash=datTimeTrace(datFile,rows,cols,rowSearch)
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
    rows=1200;
    cols=600;
end
if nargin<4
    rowSearch=cols;
end



if strfind(datFile, '.dat');   
    Fid=fopen(datFile);
status=fseek(Fid,0,1);
stackSize=floor(ftell(Fid)/(2*rows*cols)-1);
status=fseek(Fid,0,-1);
else
    error('File must be a .dat binary file in uint16 form');
end


if 1
pixelValues=fread(Fid,rows*cols,'uint16',0,'l');
initialIm=reshape(pixelValues,rows,cols);



fig=imagesc(initialIm);
display('Select area to search');
roiFlash=roipoly;
delete(fig)

end


%% search for the flash
imFlash=zeros(1,stackSize);
%cell1=imFlash;
%cell2=imFlash;
progressbar(0)
for iFrame=1:stackSize
    progressbar(iFrame/stackSize)
    pixelValues=fread(Fid,rows*cols,'uint16',0,'l');
    status=fseek(Fid,2*cols*(rows)*iFrame,-1);
if 0
imFlash(iFrame)=mean(pixelValues(roiFlash));
else
    nList=(pixelValues(roiFlash))-90;
    newVal=mean(nList(nList>median(nList)));
    if isempty(newVal);newVal=0;end
imFlash(iFrame)=newVal;
end
  %  status=fseek(Fid,2*frameNumber*1024^2,-1);

end

%%
fclose(Fid);
