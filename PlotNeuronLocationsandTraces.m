  rows=1200;cols=600;nPix=rows*cols;
  
  
  
  %% recover alignments
 alignments= load([dataFolder filesep 'alignments.mat']);
 alignments=alignments.alignments;
 
lowResFluor2BF=alignments.lowResFluor2BF;
S2AHiRes=alignments.S2AHiRes;
Hi2LowResF=alignments.Hi2LowResF;
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;



neuronList=inFrame;

%%
        timeList=[10 100 200 400];
         Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);
zOffset=offset;
rejectList=find(rejects);
%%
for iStack=timeList
    %%
 

currentFiducials=fiducialPoints{iStack};
    
    currentFiducials(rejects,:)=[];
    missingPoints=cellfun(@(x) isempty(x),currentFiducials(:,4));
    currentFiducials(missingPoints,:)={nan};
    currentFiducials=cell2mat(currentFiducials);
selectPoints=currentFiducials(neuronList,:);
selectPoints=(selectPoints);
  %  hiResIdx=find(hiResData.stackIdx==iStack)+ zOffset;
    hiResIdx=round(mean(selectPoints(:,4)));    
    inFrame=(abs(currentFiducials(:,4)-hiResIdx)<2);

    status=fseek(Fid,2*(hiResIdx(1))*nPix,-1);
    pixelValues=fread(Fid,nPix*(length(hiResIdx)),'uint16',0,'l');
    hiResImage=reshape(pixelValues,rows,cols,length(hiResIdx));
    
    segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
    activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
    segmentChannel=pedistalSubtract(segmentChannel);
    activityChannel=pedistalSubtract(activityChannel);

    figure
    imagesc(segmentChannel);colormap hot;axis square off tight
    hold on
    scatter(selectPoints(:,1),selectPoints(:,2),'+');
    
    text(selectPoints(:,1),selectPoints(:,2),[cellstr(num2str(cgIdxRev(neuronList)))],'VerticalAlignment'...
        ,'bottom', 'HorizontalAlignment','right','color',[1 1 1],...
        'fontsize',10);
    axis([100 500 100 500]);
    hold off
    savefig([dataFolder filesep 'ExampleImageRed' num2str(iStack)]);
       figure
    imagesc(activityChannel);colormap black2green;axis square off tight
    hold on
    scatter(selectPoints(:,1),selectPoints(:,2),'+');
  text(selectPoints(:,1),selectPoints(:,2),[cellstr(num2str(cgIdxRev(neuronList)))],'VerticalAlignment'...
        ,'bottom', 'HorizontalAlignment','right','color',[1 1 1],...
        'fontsize',10);
    axis([100 500 100 500]);

    hold off
       savefig([dataFolder filesep 'ExampleImageGreen' num2str(iStack)]);

    
end
%%
select=neuronList;
plotOffset=1;
topAll=0;
for class=1
    space=spaceVec(class);
   % plotOffset=i+2+plotOffset;
    plotOffset=-.5;

    
    if ~isempty(select)
%select=select(1:min(6,length(select)));

for i=1:length(select)
        plot(frameTimeTrack(FplotSelect),medfilt1(Ratio2(select(i),FplotSelect),5)'+space*(i+plotOffset),'black','linew',1)
hold on
topPoint=max(medfilt1(Ratio2(select(i),FplotSelect),5)'+space*(i+plotOffset));
topAll=max(topAll,topPoint);
%plot(R2(possibleB(i),:)'+space*i,'r')

end

%ylim([0 ethoHeight])
set(gca,'YTick',space,'fontsize',12)
xlabel('Time(s)')
ylabel('\Delta F /F0')
  text( max(frameTimeTrack(FplotSelect)+5)*ones(size(select)), space*plotOffset+[space:space:space*(i)],cellstr(num2str(cgIdxRev(select))),'VerticalAlignment'...
        ,'middle', 'HorizontalAlignment','left','color',[0 0 0],...
        'fontsize',13);
    xlim([0 max(frameTimeTrack(FplotSelect)+15)])
    %axis square
    

    end
end

       savefig([dataFolder filesep 'ExampleImageTraces']);

%%
plotOffset

