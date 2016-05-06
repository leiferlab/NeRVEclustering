function autoCreateAlignment()
%createAlignment creates an alignment mapping out of images to be used with
%worm segmentation and gcamp signal, only uses projective mapping
 

% search params

searchRad=3;
param.dim=2;
param.excessive=4;
param.quiet=1;
param.difficult=7.e2;
param.good=2;

display(['Normal convention is:'  char(10) ...
 'Red channel 2 green channel,' char(10)...
 'low mag behavior to low mag fluor,' char(10)...
 'and hiRes Red to low Res fluor'])

choice = menu('Image Setup','Split','Multiple');


%pick files
display('Select image files')
mostRecent=getappdata(0,'mostRecent');
[fileName]=uipickfiles('FilterSpec',mostRecent);
mostRecent=fileparts(fileName{1});
setappdata(0,'mostRecent',mostRecent);
display('Select initial transform');
[initialAlignment]=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration\');
if isnumeric(initialAlignment)
    initialFlag=0;
else
    initialFlag=1;
    initialAlignment=initialAlignment{1};
    initialAlignment=load(initialAlignment);
    t_concord=initialAlignment.t_concord;
end


%name alignmentfiles
%if nargin==0;
alignmentName = inputdlg('Name the alignment file:', 's');
segmentPts=initialAlignment.Sall;
activityPts=initialAlignment.Aall;
% else
%     alignment=load(alignmentName);
%     alignmentName = inputdlg('Name the alignment file:', 's');
% 
%     if choice==2
% 
%     end
% segmentPts=alignment.Sall;
% activityPts=alignment.Aall;
% end
%%
       

if choice==1
        initialIm=double(imread([fileName{1}],'tif'));

    if ~initialFlag
%% Draw 2 rectangles for the shape and activity channels
fig=imagesc(initialIm);
display('Get segmenting ROI')
rect1=getrect(gcf);
rect1=round(rect1);

display('Get Activity ROI');
rect2=getrect(gcf);
rect2=round(rect2);
    else
rect1=initialAlignment.rect1;
rect2=initialAlignment.rect2;
    end

rect1=round(rect1 +[0,0 rect1(1:2)]);
rect1(rect1<1)=1;
rect1(3)=min(rect1(3),size(initialIm,2));
rect1(4)=min(rect1(4),size(initialIm,1));

rect2=round(rect2 +[0,0 rect2(1:2)]);
rect2(rect2<1)=1;
rect2(3)=min(rect2(3),size(initialIm,2));
rect2(4)=min(rect2(4),size(initialIm,1));



nImage=length(fileName);
else

    nImage=length(fileName)/2;
end

%%
for iImage=1:nImage
    if iImage==1
        searchRad=50;
    else
        searchRad=10;
    end
    
    if choice==1

        
    initialIm=double(imread([fileName{iImage}],'tif'));

channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));

    else
    channelSegment=double(imread([fileName{iImage}],'tif'));
    channelActivity=double(imread([fileName{nImage+iImage}],'tif'));
    end
    channelSegment=bpass(channelSegment(:,:,1),1,20);
    channelActivity=bpass(channelActivity(:,:,1),1,20);
channelSegment=normalizeRange(double(channelSegment));
channelActivity=normalizeRange(double(channelActivity));

% %%
% h=imagesc(channelSegment);
% colormap(gray)
% imCon=imcontrast;
% uiwait(imCon);
% cRange=caxis;
% channelSegment(channelSegment<cRange(1))=cRange(1);
% channelSegment(channelSegment>cRange(2))=cRange(2);
% channelSegment=normalizeRange(channelSegment);
% 
% %%
% h=imagesc(channelActivity);
% imCon=imcontrast;
% uiwait(imCon);
% colormap(gray)
% cRange=caxis;
% channelActivity(channelActivity<cRange(1))=cRange(1);
% channelActivity(channelActivity>cRange(2))=cRange(2);
% channelActivity=normalizeRange(channelActivity);
% close all


%%
activityBW=channelActivity>graythresh(channelActivity);
segmentBW=channelSegment>graythresh(channelSegment);
activityBW=AreaFilter(activityBW,20);
segmentBW=AreaFilter(segmentBW,20);
statsActivity=regionprops(activityBW,channelActivity,'WeightedCentroid');
statsSegment=regionprops(segmentBW,channelSegment,'WeightedCentroid');

aCentroids=cell2mat({statsActivity.WeightedCentroid}');
sCentroids=cell2mat({statsSegment.WeightedCentroid}');

if initialFlag
    aCentroids2=transformPointsForward(t_concord,...
        aCentroids);
    sCentroids2=sCentroids;
    if iImage==1
    aCentroids2=bsxfun(@minus,aCentroids2,mean(aCentroids2,1));
    sCentroids2=bsxfun(@minus,sCentroids,mean(sCentroids,1));

    aCentroids2=aCentroids2+1000;
    sCentroids2=sCentroids2+1000;
    end
    
        
else
    [activityT,segmentT]=cpselect(channelActivity,channelSegment,...
                activityPts,segmentPts,'Wait',true);
    t_concord = fitgeotrans(activityT,segmentT,'projective');
    aCentroids2=transformPointsForward(t_concord,...
        aCentroids);
    sCentroids2=sCentroids;
    initialFlag=1;
end


trackInput=[aCentroids2 (1:size(aCentroids2,1))' ones(size(aCentroids2(:,1))) ; ...
    sCentroids2 (1:size(sCentroids2,1))' 2*ones(size(sCentroids2(:,1)))];
trackInput=trackInput(all(trackInput>0,2),:);
TrackOut=nan;
for iSearch=0:searchRad
    if all(isnan(TrackOut(:)))
        TrackOut=trackJN(trackInput,searchRad-iSearch,param);
    end
end
aTrack=TrackOut(TrackOut(:,4)==1,3);
sTrack=TrackOut(TrackOut(:,4)==2,3);

activityPts=cat(1,activityPts,aCentroids(aTrack,:));
segmentPts=cat(1,segmentPts,sCentroids(sTrack,:));

t_concord = fitgeotrans(activityPts,segmentPts,'projective');

activityPts2=transformPointsForward(t_concord,...
activityPts);

scatter(activityPts2(:,1),activityPts2(:,2),'x')
hold on
scatter(segmentPts(:,1),segmentPts(:,2))
plot([segmentPts(:,1) activityPts2(:,1)]',...
    [segmentPts(:,2) activityPts2(:,2)]');
hold off
pause(.1)
%% Select control points and create transform

            
end
%%
%reject points taht are more than 2 sd further than average
zDistance=zscore(sqrt(sum((activityPts2-segmentPts).^2,2)));
zSelect=zDistance<2;
Aall=activityPts(zSelect,:);
Sall=segmentPts(zSelect,:);
t_concord = fitgeotrans(Aall,Sall,'projective');
Rsegment = imref2d(size(channelSegment));
activityRegistered = imwarp(channelActivity,t_concord,'OutputView',Rsegment);
padRegion=activityRegistered==0;
padRegion=imdilate(padRegion,true(3));
%% save 
if choice==1
save(['Y:\CommunalCode\3dbrain\registration\' alignmentName{1}],'rect1','rect2','t_concord'...
    ,'Rsegment','padRegion','initialIm','Sall','Aall','fileName')
else
save(['Y:\CommunalCode\3dbrain\registration\' alignmentName{1}],'t_concord'...
    ,'Rsegment','Sall','Aall','fileName')
    
end





