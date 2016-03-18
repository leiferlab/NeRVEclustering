%%% flatCheck takes a set of .dat file movies from the sCMOS camera and
%%% calculates the radial inhomogeneity of the filed of view by taking
%%% circular projections of the average image.

 % select folders with image files
 display('select folders with image files')
folderList=uipickfiles('FilterSpec','O:');
imsize=[1200 600];
method='mean';
%% loop through files and take mean projections
for iFolder=1:length(folderList); 
    bfFile=[folderList{iFolder} filesep 'sCMOS_Frames_U16_1024x1024.dat'];
imageProj=zProjectDatFile(bfFile,imsize,method,3); %take mean projection
projAll{iFolder}=imageProj;
end

%% load alignment files
    alignments=load([dataFolder filesep 'alignments']);
    alignments=alignments.alignments;
lowResFluor2BF=alignments.lowResFluor2BF;
S2AHiRes=alignments.S2AHiRes;
Hi2LowResF=alignments.Hi2LowResF;
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;

%% create mask seperate two halves of the DV2 image


cropRegion=rect1;
cropPointsx=[cropRegion(1),cropRegion(1),cropRegion(1)+cropRegion(3),...
    cropRegion(1)+cropRegion(3)];
cropPointsy=[cropRegion(2),cropRegion(2)+cropRegion(4),...
    cropRegion(2)+cropRegion(4),cropRegion(2)];

cropMask=poly2mask(cropPointsx,cropPointsy,imsize(1),imsize(2));
imMeans=cellfun(@(x) mean(x(cropMask)),projAll);
imSTD=cellfun(@(x) std(x(cropMask)),projAll);
cropProj=imcrop(projAll{1},cropRegion);

%% plot circular projections for each movie
close all
for iFolder=1:length(folderList)
cropProj=imcrop(projAll{iFolder},cropRegion); %crop image
%create cartesian coordinate system
[cartX,cartY]=meshgrid(1:size(cropProj,1),1:size(cropProj,2));
cartX=cartX-mean(cartX(:));
cartY=cartY-mean(cartY(:));
%creat polar coordinates
[phi,r]=cart2pol(cartX,cartY);
r=round(r+.5);
% radial project with accumarray
rplot=accumarray(r(:),cropProj(:),[],@mean)./imMeans(iFolder);
plot(rplot);
hold on
end

