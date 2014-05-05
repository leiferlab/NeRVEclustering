

%load folder and extract and align timing data.
imFolder=uigetdir;
imNames=dir([imFolder filesep '*.tif']);
alignFlag=0;
if alignFlag
%% better initial Im by averaging multiple stacks and doing a max projection
%!!! may be better to just take images from BEFORE the flash,
%photobleaching appears to make these images better than my averaging. 
clear worm

for iImage=1

        temp=double(imread([imFolder filesep imNames(iImage).name],'tif'));
                temp=pixelIntensityCorrection(temp);
                if iImage==1
                    initialIm=temp;
                else
                    
        initialIm=initialIm+(temp);
                end
end





%% Draw 2 rectangles for the shape and activity channels
fig=imagesc(initialIm);
display('Get segmenting ROI')
rect1=getrect(gcf);
rect1=round(rect1);
rectSize1=rect1(3:4);
rect1=round(rect1 +[0,0 rect1(1:2)]);
channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
display('Get Activity ROI');
rect2=getrect(gcf);
rect2=round(rect2);
rectSize2=rect2(3:4);
rect2=round(rect2 +[0,0 rect2(1:2)]);
channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));


channelSegment=normalizeRange(double(channelSegment));
channelActivity=normalizeRange(double(channelActivity));
close all
%% Select control points and create transform
[activityPts,segmentPts]=cpselect(wiener2(channelActivity,[3,3],1),channelSegment,...
                'Wait',true);
t_concord = fitgeotrans(activityPts,segmentPts,'projective');
Rsegment = imref2d(size(channelSegment));
activityRegistered = imwarp(channelActivity,t_concord,'OutputView',Rsegment);
padRegion=activityRegistered==0;
padRegion=imdilate(padRegion,true(3));
%% save 
save([imFolder filesep 'stackInfo'],'rect1','rect2','t_concord'...
    ,'Rsegment','rectSize1','rectSize2','padRegion','imFolder','initialIm');

else
    [regFile,regFolder]=uigetfile('Y:\CommunalCode\3dbrain\registration\');
    load([regFolder filesep regFile]);
end

if alignFlag==2
    
    save(['Y:\CommunalCode\3dbrain\registration\' datestr(date,29)],'rect1','rect2','t_concord'...
    ,'Rsegment','rectSize1','rectSize2','padRegion','initialIm')
end


%%

%% segment subimages and create masks
[folderList] = uipickfiles;
for iFolder=1:length(folderList)
imFolder=folderList{iFolder};
imNames=dir([imFolder filesep '*.tif']);
mkdir([imFolder filesep 'stackData']);

movieLength=length(imNames);
progressbar(0)
for iImage=1:movieLength;
    progressbar(iImage/movieLength);
    tic

    % load stacks

        temp=double(imread([imFolder filesep imNames(iImage).name],'tif'));
        temp=pixelIntensityCorrection(temp);
        temp_activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
        worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
        temp_activity=imwarp(temp_activity,t_concord,'OutputView',Rsegment);
        temp_activity(padRegion)=median(temp_activity(~padRegion));
    %    activity=bpass_jn((temp_activity),1,[20,20]);
        activity=temp_activity;
    imsize=size(worm);  
  %do segmentation
    wormMask=WormSegmentHessian2D(worm);
   wormMask= bwmorph(wormMask,'clean');
    %look up intensities on both channels, after a bit of dilation
    wormLabelMask=imdilate(bwlabeln(wormMask),true(5,5));
wormcc=bwconncomp(wormMask);
stats=regionprops(wormcc,'Centroid','Area');
centroids=reshape([stats.Centroid],2,[])';
Rintensities=cellfun(@(x) trimmean(worm(x),20),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) trimmean(activity(x),20),[wormcc.PixelIdxList])';

    %interpolate Z properly and scale
Volume=[stats.Area]';
%save outputs in unique file
outputFile=[imFolder filesep 'stackData' filesep 'stack' num2str(iImage,'%04d') 'data'];

save(outputFile,'centroids','Rintensities','Gintensities','Volume',...
    'wormMask');
display(['Completed stack' num2str(iImage,'%04d') ' in ' num2str(toc) ' seconds']);
end
end
