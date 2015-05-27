%dataFolder='E:\20141212\BrainScanner20141212_145951\';
dataFolder=uipickfiles;
dataFolder=dataFolder{1};

imSize=[1200 600];
%[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.1 0.01], [0.1 0.01]);





%%
imageFolder=[dataFolder filesep 'rawVideoFeeds'];

mkdir(imageFolder);
%%

%%
%hiResRegistration=load('Y:\CommunalCode\3dbrain\registration\20141212HiResS2A.mat');
hiResRegistration=uipickfiles('filterspec','Y:\CommunalCode\3dbrain\registration\');
 hiResRegistration=load(hiResRegistration{1});
rect1=hiResRegistration.rect1;
    rect2=hiResRegistration.rect2;
    t_concord=hiResRegistration.t_concord;
    Rsegment=hiResRegistration.Rsegment;
    padRegion=hiResRegistration.padRegion;
    
%%

Pmap=cbrewer('seq','Blues',100);
Bmap=cbrewer('seq','YlOrRd',100);
Fmap=cbrewer('seq','YlGnBu',100);
Fmap=flipud(Fmap);
%%
fclose all;

    Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat'] );
%frameRange=find(hiResData.stackIdx==stackList(1),1,'first'):find(hiResData.stackIdx==stackList(end),1,'last');
frameRange=1000:5000;

;%%
pointHistory=[];
frameCounter=1;



for iStack=1:length(frameRange)
    tic
  
        %%
        try
    hiResIdx=frameRange(iStack);
    imName=['image' num2str(iStack,'%3.5d') '.tif'];
    matName=['image' num2str(iStack,'%3.5d') '.mat'];
%     metaData=load([metaFolder filesep matName]);
%     metaData=metaData.metaData.metaData;
 %   wormInfo=imfinfo([hiResSegmentFolder filesep imName]);
%     
% worm=stackLoad([hiResSegmentFolder filesep imName]);
% activity=stackLoad([hiResActivityFolder filesep  imName]);
% if metaData.zVoltage(1)>metaData.zVoltage(end);
%     worm=flipdim(worm,3);
%     activity=flipdim(activity,3);
% end
% worm=normalizeRange(worm);
% activity=normalizeRange(activity);


    

    status=fseek(Fid,2*hiResIdx*imSize(1)*imSize(2),-1);
    temp=fread(Fid,imSize(1)*imSize(2),'uint16',0,'l');
    temp=reshape(temp,imSize);
    worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
    
        activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
        activity=imwarp(activity,t_concord,'OutputView',Rsegment,'fillvalue',nan);
   worm=pedistalSubtract(worm);
        activity=pedistalSubtract(activity);
    worm(worm==0)=nan;

    activity(activity==0)=nan;

%%
    filename=[imageFolder filesep 'frame' num2str(iStack,'%3.5d')];
 %   filename2=[imageFolder2 filesep 'frame' num2str(frameCounter,'%3.5d')];
   % worm
    subplot(1,2,1)
    
    imagesc(worm);
    colormap(hot)
colorbar off
    set(gcf,'color','black');
    caxis([0 200]);
    freezeColors
axis equal off tight

%eigenPosture
     subplot(1,2,2);
% %delete(gca)
% colorIdx=interp1(colorLookup,1:256,bfeigenZ(bfIdx(iSlice)),'nearest');
%     colorVector=colorMap(colorIdx,:);
%     scatter3(bfeigenProj(bfIdx(iSlice),1),bfeigenProj(bfIdx(iSlice),2),...
%         bfeigenZ(bfIdx(iSlice)),'markerFacecolor',color,'markeredgecolor','none');        
%    xlim([-8,8]);
%     ylim([-8,8]);
%     zlim([-20 20]);
%     xlabel('Eig1');ylabel('Eig2');zlabel('Eig3');
%     
       imagesc(activity)
    colormap(pmkmp(123,'LinearL'));
        caxis([0 200]);

colorbar off
   % imagesc(worm(:,:,iSlice));
    %set(gcf,'color','w');
 %   caxis([0 1]);
    freezeColors
axis equal off tight
%caxis([0,1]);


    
    
  %  xlabel('x (mm)');
%     hold on
% pointHistory=[pointHistory;[xPosFrame(bfIdx(iSlice)),xPosFrame(bfIdx(iSlice))]];
%  historyLength=size(pointHistory,1);
% plot(pointHistory(1:10:end,1)...
%     ,pointHistory(1:10:end,2))
% hold off


set(gcf,'color',[0 0 0])
drawnow
saveas(gcf,filename,'jpeg');



    display(['completed stack:' num2str(iStack) ' in ' num2str(toc) 's']);
        catch ME
      display(['error stack:' num2str(iStack) ' in ' num2str(toc) 's']);
%save(['errorStack' num2str(iStack)], 'ME') 
        end
        
end
