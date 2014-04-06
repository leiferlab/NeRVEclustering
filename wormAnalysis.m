
imFolder=uigetdir;
[initialIm,stackInfo]=timeTraceAnalysis(imFolder);

%%
fig=imagesc(initialIm);
rect1=getrect(gcf);
rect1=round(rect1);
rectSize1=rect1(3:4);
rect1=round(rect1 +[0,0 rect1(1:2)]);
channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

rect2=getrect(gcf);
rect2=round(rect2);
rectSize2=rect2(3:4);
rect2=round(rect2 +[0,0 rect2(1:2)]);
channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));


channelSegment=normalizeRange(double(channelSegment));
channelActivity=normalizeRange(double(channelActivity));

%%
[activityPts,segmentPts]=cpselect(channelActivity,channelSegment,...
                'Wait',true);
t_concord = fitgeotrans(activityPts,segmentPts,'projective');
Rsegment = imref2d(size(channelSegment));
activityRegistered = imwarp(channelActivity,t_concord,'OutputView',Rsegment);
            
%%
for i=2
stackSize=length(stackInfo(i).fileNames);
    worm=zeros(rectSize1(2),rectSize1(1),stackSize);
    activity=worm;
    for slice=1:stackSize
        
        temp=imread([imFolder filesep stackInfo(i).fileNames(slice).name],'tif');
        temp_activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
        worm(:,:,slice)=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
        activity(:,:,slice)=imwarp(temp_activity,t_concord,'OutputView',Rsegment);
    end
    imsize=size(worm);  
    worm=image_resize(worm,imsize(1),imsize(2),2*imsize(3));
    
end
