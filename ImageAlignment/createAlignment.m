function createAlignment(alignmentName)
%createAlignment creates an alignment mapping out of images to be used with
%worm segmentation and gcamp signal, only uses projective mapping
 
display(['Normal convention is:'  char(10) ...
 'Hi Res Red channel 2 green channel,' char(10)...
 'low mag behavior to low mag fluor,' char(10)...
 'and hiRes Red to low Res fluor'])

choice = menu('Image Setup','Align 2 halves','Align different files');


%pick files
[fileName]=uipickfiles('FilterSpec','E:\');
%name alignmentfiles
if nargin==0;
alignmentName = inputdlg(...
    'Name the alignment file, specify the date and type of alignment :', 's');
segmentPts=[];
activityPts=[];
else
    alignment=load(alignmentName);
    alignmentName = inputdlg(...
    'Name the alignment file, specify the date and type of alignment :', 's');


    if choice==2

    end
segmentPts=alignment.Sall;
activityPts=alignment.Aall;
end
%%
       

if choice==1
        initialIm=double(imread([fileName{1}],'tif'));

    if nargin==0

%% Draw 2 rectangles for the shape and activity channels
fig=imagesc(initialIm);
display('Get segmenting ROI, normally zero to middle')
rect1=getrect(gcf);
rect1=round(rect1);

display('Get Activity ROI, normally middle to end');
rect2=getrect(gcf);
rect2=round(rect2);
    else
rect1=alignment.rect1;
rect2=alignment.rect2;
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
    if choice==1

        
    initialIm=double(imread([fileName{iImage}],'tif'));

channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));

    else
    channelSegment=double(imread([fileName{iImage}],'tif'));
    channelActivity=double(imread([fileName{nImage+iImage}],'tif'));
    end
    
channelSegment=normalizeRange(double(channelSegment(:,:,1)));
channelActivity=normalizeRange(double(channelActivity(:,:,1)));

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
if ~isempty(segmentPts)
[activityPts,segmentPts]=cpselect(channelActivity,channelSegment,...
                activityPts,segmentPts,'Wait',true);
else
    [activityPts,segmentPts]=cpselect(channelActivity,channelSegment,...
                'Wait',true);
end

            
end
%%
Aall=activityPts;
Sall=segmentPts;
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



end

