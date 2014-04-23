
parent=uigetdir;
folders=dir(parent);
folders=folders([folders.isdir]);
folders=folders(3:end);
pixI=[];
pixVar=[];
for iFolder=1:length(folders)

imFolder=folders(iFolder).name;
images=dir([parent filesep imFolder filesep '*.tif']);


    imAll=double(imread([parent filesep imFolder filesep images(1).name],'tif'));
for iImage=2:length(images);
    imAll(:,:,iImage)=double(imread([parent filesep imFolder filesep images(iImage).name],'tif'));
    progressbar(iFolder/length(folders),iImage/length(images));
end

pixI(:,:,iFolder)=mean(imAll,3);
pixVar(:,:,iFolder)=var(imAll,[],3);

end 
%% do least squares fit
pixB=pixI(:,:,1);
pixI=bsxfun(@minus,pixI,min(pixI,[],3));
varpixI=var(pixI,[],3);
varpixVar=var(pixVar,[],3);
covarPix=sum(bsxfun(@minus,pixI,mean(pixI,3)).*bsxfun(@minus,pixVar,mean(pixVar,3)),3);
pixM=covarPix./varpixI;

pixM=pixM./trimmean(pixM(:),30);

save([parent filesep 'PixelCalibration'],'pixM','pixB');
