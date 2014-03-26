%% load excel file, take outputs and find peak
imFolder=uigetdir();
p=dir([imFolder filesep '*.xlsx']);
p=p(1).name;
imFiles=dir([imFolder filesep '*tif']);
output=xlsread([imFolder filesep p],1,'A:C');

%%
startWave=output(:,3);
startWave=startWave-mean(startWave);
startWave=cumsum(startWave>.2)>1;

imageWave=output(:,2)>mean(output(:,2));
zWave=output(:,1);
zWave=smooth(zWave,100);
imageRise=[0;diff(single(imageWave))>.2];
imageFall=[diff(single(imageWave))<-.2;0];
[x,lags]=xcorr(imageFall,imageRise,100);
shift=lags(x>range(x)/2 & lags'>0);
shift=shift(1);
imageSpikes=circshift(imageRise,round(shift/2));


imageZ=zWave(imageSpikes & startWave);
zgrad=medfilt1(gradient(zWave,5),5)>0;
zgrad=zgrad(imageRise & startWave);
stackIdx=[0;cumsum(abs(diff(zgrad)))];



%% load images and find flash in images using user defined ROI

stackSize=length(imFiles);
initialIm=(imread([imFolder filesep imFiles(1).name], 'tif'));
progressbar;
for i=round(stackSize/2):round(stackSize/2)+60;
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    initialIm=max(initialIm,temp);
    progressbar((i-round(stackSize/2))/60);

end
imsize=size(initialIm);
fig=imagesc(initialIm);
display('Select area to find flash');
roiFlash=roipoly;
delete(fig)

fig=imagesc(initialIm);
display('Select an cell to try to align');
roiCell1=roipoly;
delete(fig)

fig=imagesc(initialIm);
display('Select another');
roiCell2=roipoly;
close all
%% search for the flash
imFlash=zeros(1,stackSize);
cell1=imFlash;
cell2=imFlash;
progressbar;
for i=1:10:stackSize/10
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    imFlash(i)=mean(temp(roiFlash));
    cell1(i)=mean(temp(roiCell1));
    cell2(i)=mean(temp(roiCell2));
    progressbar(i/stackSize*10);
end
progressbar(1);
flashMean=mean(imFlash(imFlash>0));
flashStd=std(imFlash(imFlash>0));


while ~any(imFlash>(flashMean+5*flashStd))
    i=i+10;
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    imFlash(i)=mean(temp(roiFlash));
    flashMean=mean(imFlash(imFlash>0));
flashStd=std(imFlash(imFlash>0));
    
end
i=find(imFlash==max(imFlash));
for k=max(i-100,0):min(i+100,stackSize)
       temp=(imread([imFolder filesep imFiles(k).name], 'tif'));
    imFlash(k)=mean(temp(roiFlash)); 
        progressbar((k-i+100)/200);

    
end
imFlash(imFlash==0)=flashMean;

%%

flashWave=imFlash;
flashWave=flashWave-min(flashWave);
flashWave=flashWave>(mean(flashWave)*5);
flashWave=[0,diff(flashWave)>.5];
startFlash=find(flashWave>.5);

if length(startFlash)>1
    endFlash=startFlash(end);
    startFlash=startFlash(1);
else
    endFlash=stackSize;
end

images=imFiles(startFlash:endFlash);


%%
fig=imagesc(initialIm);
rect=getrect(gcf);
rect=round(rect);
rectSize=rect(3:4);
rect=round(rect +[0,0 rect(1:2)]);

%%
for i=2
    imSelect=find(stackIdx==i);
    hyper_stack=zeros(rectSize(2),rectSize(1),length(imSelect));
    for slice=1:length(imSelect)
        
        temp=imread([imFolder filesep images(imSelect(slice)).name],'tif');
        
        hyper_stack(:,:,slice)=temp((rect(2)+1):rect(4),(1+rect(1)):rect(3));
    end
end








    



