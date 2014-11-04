rows=1200;
cols=600;
folders=uipickfiles;
nPix=rows*cols;
%%
progressbar(0,0);
%pixI=zeros(rows,cols,length(folders));
%pixVar=zeros(rows,cols,length(folders));
    folderCounter=1;

for iFolder=1:length(folders)
    
Fid=fopen([folders{iFolder}  filesep 'sCMOS_Frames_U16_1024x1024.dat']);
if Fid~=-1
fseek(Fid,0,1);
stackSize=floor(ftell(Fid)/(2*rows*cols)-1);
fseek(Fid,0,-1);
if stackSize~=0
imAll=uint16(zeros(rows*cols,stackSize));
for iImage= 1 :stackSize 
    imAll(:,iImage)=fread(Fid,nPix,'uint16',0,'l');
    
    progressbar(iFolder/length(folders),iImage/stackSize);
end
fclose(Fid);

pixI(:,:,folderCounter)=reshape(mean(imAll,2),rows,cols);
pixVar(:,:,folderCounter)=reshape(var(double(imAll),[],2),rows,cols);
folderCounter=folderCounter+1;
end
end
end 

nanmap=squeeze(pixI(1,1,:));
pixI=pixI(:,:,~isnan(nanmap));
pixVar=pixVar(:,:,~isnan(nanmap));

meanmap=squeeze(nanmean(nanmean(pixI,1),2));
[~,ib]=sort(meanmap);
pixI=pixI(:,:,ib);
pixVar=pixVar(:,:,ib);


%% do least squares fit


pixB=pixI(:,:,1);
pixI=bsxfun(@minus,pixI,min(pixI,[],3));
varpixI=var(pixI,[],3);
varpixVar=var(pixVar,[],3);
covarPix=sum(bsxfun(@minus,pixI,mean(pixI,3)).*bsxfun(@minus,pixVar,mean(pixVar,3)),3);
pixM=covarPix./varpixI;
pixM(pixM<0)=0;
pixM=pixM./trimmean(pixM(:),30);

save([fileparts(folders{1}) filesep 'PixelCalibration'],'pixM','pixB');
