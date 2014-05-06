function createAlignment()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[fileName,pathName]=uigetfile('*.tif','MultiSelect','on');




    initialIm=double(imread([pathName filesep fileName{1}],'tif'));

%% Draw 2 rectangles for the shape and activity channels
fig=imagesc(initialIm);
display('Get segmenting ROI')
rect1=getrect(gcf);
rect1=round(rect1);
rectSize1=rect1(3:4);
rect1=round(rect1 +[0,0 rect1(1:2)]);

display('Get Activity ROI');
rect2=getrect(gcf);
rect2=round(rect2);
rectSize2=rect2(3:4);
rect2=round(rect2 +[0,0 rect2(1:2)]);

Aall=[];
Sall=[];
for iImage=1:length(fileName)
    initialIm=double(imread([pathName filesep fileName{iImage}],'tif'));

channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));


channelSegment=normalizeRange(double(channelSegment));
channelActivity=normalizeRange(double(channelActivity));
close all
%% Select control points and create transform
if ~isempty(Aall)
[activityPts,segmentPts]=cpselect(channelActivity,channelSegment,...
                Aall,Sall,'Wait',true);
else
    [activityPts,segmentPts]=cpselect(channelActivity,channelSegment,...
                'Wait',true);
end
 Aall=cat(1,Aall,activityPts);
 Sall=cat(1,Sall,segmentPts);
            
end


t_concord = fitgeotrans(Aall,Sall,'projective');
Rsegment = imref2d(size(channelSegment));
activityRegistered = imwarp(channelActivity,t_concord,'OutputView',Rsegment);
padRegion=activityRegistered==0;
padRegion=imdilate(padRegion,true(3));
%% save 
fileDateString=fileparts(fileparts(pathName));
fileDateString=fileDateString(strfind(fileDateString,filesep)+1:end);
save(['Y:\CommunalCode\3dbrain\registration\' fileDateString],'rect1','rect2','t_concord'...
    ,'Rsegment','rectSize1','rectSize2','padRegion','initialIm')



end

