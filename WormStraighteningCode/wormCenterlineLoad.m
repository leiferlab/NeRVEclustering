%%

smooth_weight = .001; 
image_weight = 100; 
delta_t = 4; 
margin = 10; 



%% select file name

movieFile=uipickfiles('FilterSpec','Y:\PanNeuronal\20140820');
movieFile=movieFile{1};
centerLineFile=uipickfiles('FilterSpec','Y:\PanNeuronal\20140820');

%% load vidObj using VideoReader
if strfind(movieFile,'.avi')
vidObj = VideoReader(movieFile);
lastFrame = read(vidObj, inf);
 numFrames= vidObj.NumberOfFrames;
 aviFlag=1;
elseif isdir(movieFile)
    movFiles=dir([movieFile filesep '*.jpeg']);
    aviFlag=0;
    
end

    
    
    
    
 load(centerLineFile{1})

 
 %%
 stretchSize=15;
 figure
 

 
for iFrame=1:size(centerline,3);
    if aviFlag
    lastFrame = read(vidObj,iFrame);
lastFrame=normalizeRange(sum(double(lastFrame),3));
    else
lastFrame=imread([movieFile filesep movFiles(iFrame).name]);
lastFrame=normalizeRange(sum(double(lastFrame),3));

    end
CL=centerline(:,:,iFrame);
CL=[interp1(CL(:,1),-stretchSize+1:100+stretchSize,'*linear','extrap')',...
    interp1(CL(:,2),-stretchSize+1:100+stretchSize,'*linear','extrap')'];
    [newIm,newX,newY]=wormStraightening(CL,lastFrame,60);

    
    
%     imagesc(sum(lastFrame,3))
%     hold on
%     
%     plot(centerline(:,2,iFrame),centerline(:,1,iFrame))
%     drawnow
%     hold off
subplot(1,2,1);
imagesc(lastFrame);
hold on
 plot(CL(:,2),CL(:,1))
  plot(newY',newX','black')
  hold off

subplot(1,2,2);

imagesc(newIm)
%     if iFrame==1
% phi = zeros(size(newIm)); 
% phi(margin:end-margin, margin:end-margin) = 1; 
% phi = ac_reinit(phi-.5); 
%     end
%     phi = ac_ChanVese_model(max(phi2(:))*normalizeRange(newIm), phi, smooth_weight, image_weight, delta_t, 1);
%     
% hold on
%   [C,h]=contour(phi2,[0,0],'w');axis equal;
%   hold off
%   
  
    axis equal
    drawnow
end

%%