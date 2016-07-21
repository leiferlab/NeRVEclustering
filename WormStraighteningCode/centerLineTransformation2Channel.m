%%

smooth_weight = .001; 
image_weight = 100; 
delta_t = 4; 
margin = 10; 
%% load registration
display('Select Registration mat File');
[regFile,regFolder]=uigetfile('Y:\CommunalCode\3dbrain\registration\');
    load([regFolder filesep regFile]);

%% select file name

behaviorMovie=uipickfiles('FilterSpec','Y:\PanNeuronal\20140905\PNWorm1\');
behaviorMovie=behaviorMovie{1};
fluorMovie=uipickfiles('FilterSpec','Y:\PanNeuronal\20140905\PNWorm1\');
fluorMovie=fluorMovie{1};
%%
centerLineFile=uipickfiles('FilterSpec','Y:\PanNeuronal\20140905\PNWorm1\');
fluor2bfIdx=YamlFlashAlign;
 load(centerLineFile{1})

%% load vidObj using VideoReader
if strfind(behaviorMovie,'.avi')
behaviorVidObj = VideoReader(behaviorMovie);
lastFrame = read(behaviorVidObj, inf);
 numFrames= behaviorVidObj.NumberOfFrames;
 
 fluorVidObj= VideoReader(fluorMovie);
 
 
 aviFlag=1;
elseif isdir(behaviorMovie)
    movFiles=dir([behaviorMovie filesep '*.jpeg']);
    aviFlag=0;
    
end
 
 %%
 stretchSize=15;
 figure
 

 firstFrame=find(~isnan(fluor2bfIdx),1,'first');
for iFrame=firstFrame:size(centerline,3);
    iFrameFluor=fluor2bfIdx(iFrame);
    if aviFlag
    bfFrame = read(behaviorVidObj,iFrame);
bfFrame=normalizeRange(sum(double(bfFrame),3));
fluorFrame=read(fluorVidObj,iFrameFluor);
fluorFrame=normalizeRange(sum(double(fluorFrame),3));

        fluorFrame=imwarp(fluorFrame,t_concord,'OutputView',Rsegment);


    else
bfFrame=imread([behaviorMovie filesep movFiles(iFrame).name]);
bfFrame=normalizeRange(sum(double(bfFrame),3));

    end
CL=centerline(:,:,iFrame);
CL=[interp1(CL(:,1),-stretchSize+1:100+stretchSize,'*linear','extrap')',...
    interp1(CL(:,2),-stretchSize+1:100+stretchSize,'*linear','extrap')'];
CL2=bsxfun(@plus,CL,[0,80]);


    [newBF,newX,newY]=wormStraightening(CL,bfFrame,60);
    [newFluor,~,~]=wormStraightening(CL2,fluorFrame,60);

    
    
%     imagesc(sum(lastFrame,3))
%     hold on
%     
%     plot(centerline(:,2,iFrame),centerline(:,1,iFrame))
%     drawnow
%     hold off
subplot(2,2,[1]);
imagesc(fluorFrame);
hold on
 plot(CL2(:,2),CL2(:,1))
 % plot(newY',newX','black')
  hold off
axis equal
  subplot(2,2,3)
  imagesc(bfFrame);
hold on
 plot(CL(:,2),CL(:,1))
  plot(newY',newX','black')
  hold off
axis equal
subplot(2,2,4);

imagesc(newBF);
axis equal
subplot(2,2,2);
imagesc(newFluor);
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