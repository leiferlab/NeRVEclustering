function createAlignment()
%createAlignment creates an alignment mapping out of images to be used with
%worm segmentation and gcamp signal

choice = menu('Image Setup','Split','Multiple');



%pick files
[fileName]=uipickfiles();
%name alignmentfiles
alignmentName = inputdlg('Name the alignment file:', 's');


if choice==1
    initialIm=double(imread([fileName{1}],'tif'));

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
rect1(rect1<1)=1;
rect2(rect2<1)=1;

nImage=length(fileName);
else
    channelSegment=double(imread([fileName{1}],'tif'));
    channelActivity=double(imread([fileName{2}],'tif'));
    nImage=1;
end


Aall=[];
Sall=[];
for iImage=1:nImage
    if choice==1

        
    initialIm=double(imread([fileName{iImage}],'tif'));

channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));

    end
channelSegment=normalizeRange(double(channelSegment));
channelActivity=normalizeRange(double(channelActivity));

%%
h=imagesc(channelSegment);
colormap(gray)
imCon=imcontrast;
uiwait(imCon);
cRange=caxis;
channelSegment(channelSegment<cRange(1))=cRange(1);
channelSegment(channelSegment>cRange(2))=cRange(2);
channelSegment=normalizeRange(channelSegment);

%%
h=imagesc(channelActivity);
imCon=imcontrast;
uiwait(imCon);
colormap(gray)
cRange=caxis;
channelActivity(channelActivity<cRange(1))=cRange(1);
channelActivity(channelActivity>cRange(2))=cRange(2);
channelActivity=normalizeRange(channelActivity);
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
%%

t_concord = fitgeotrans(Aall,Sall,'projective');
Rsegment = imref2d(size(channelSegment));
activityRegistered = imwarp(channelActivity,t_concord,'OutputView',Rsegment);
padRegion=activityRegistered==0;
padRegion=imdilate(padRegion,true(3));
%% save 
if choice==1
save(['Y:\CommunalCode\3dbrain\registration\' alignmentName{1}],'rect1','rect2','t_concord'...
    ,'Rsegment','rectSize1','rectSize2','padRegion','initialIm')
else
save(['Y:\CommunalCode\3dbrain\registration\' alignmentName{1}],'t_concord'...
    ,'Rsegment','padRegion')
    
end



end

