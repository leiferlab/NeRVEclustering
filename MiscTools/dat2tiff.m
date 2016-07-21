
%loads dat file and calculates the average intensity of each pixel. method
%is not yet available. image must be [1200 600], and be uint16. 
bgFile=uipickfiles;
    imsize=[1200 600];
skip=1;
nPix=prod(imsize);
bgFile=bgFile{1};
dataFolder=fileparts(bgFile);
Fid=fopen(bgFile);

status=fseek(Fid,0,1);
lastFrame=ftell(Fid)/(prod(imsize)*2);

%%
status=fseek(Fid,0,-1);
camFrameData=importdata([dataFolder filesep 'CameraFrameData.txt']);
camFrameData=camFrameData.data;
frameNum=camFrameData(:,2);
saveFlag=find(diff(frameNum));
saveBreaks=[(diff([saveFlag])>1);1];
saveLength=diff([0;find(saveBreaks)]);
imageName=inputdlg('Image Names');
for iChunk=1:length(saveLength);
    pixelValues=fread(Fid,prod(imsize)*saveLength(iChunk),'uint16',0,'l');
Im=mean(reshape(pixelValues,imsize(1),imsize(2),saveLength(iChunk)),3);
outputImage=[dataFolder filesep imageName num2str(iChunk,'%3.2d') '.tif'];
tiffwrite(outputImage,Im)

end

