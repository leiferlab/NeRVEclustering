
%%

maxFrame=43326; %%% CHECK FOR LAST USEFUL FRAME IN VIDEO


cline_para.tipRegion=20;
cline_para.endRepulsion=.3;
cline_para.repulsionD=20;
cline_para.heat=1;
cline_para.CLalpha=5;
cline_para.CLbeta=70;
cline_para.gamma=10;
cline_para.kappa=15;
cline_para.endkappa=20;
cline_para.gradient_force=50;
cline_para.showFlag=0;
cline_para.iterations=400;
cline_para.stretching_force_factor=[2 2];
cline_para.refSpring=1;
cline_para.stretch_ends_flag=1;
cline_para.refSpring=2;
cline_para.springForce=1.01;
cline_para.refL=6;
cline_para.memForce=.02;
close all
gaussFilter=fspecial('gaussian',20,75);
%%
%% set up high mag videos
rows=1200;
cols=600;
nPix=rows*cols;

[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);



%% set up different kernals

gaussKernal2=gausswin(200);
gaussKernal2=convnfft(gaussKernal2,gaussKernal2');



Sfilter=max(gaussKernal2(:))-gaussKernal2;

Sfilter(Sfilter<.1)=-(.1-Sfilter(Sfilter<.1))*80;%Sfilter(Sfilter<.01)-.3;

Sfilter(Sfilter>.8)=0;
Sfilter(Sfilter>0)=1;%nnz(Sfilter<0)/nnz(Sfilter>0);
%Sfilter=Sfilter-mean(Sfilter(:));

Sfilter2=max(gaussKernal2(:))-gaussKernal2;
Sfilter2(Sfilter2<.01)=Sfilter2(Sfilter2<.01)-.3;
Sfilter2(Sfilter2>.6)=0;

Sfilter2(Sfilter2>0)=nnz(Sfilter2<0)/nnz(Sfilter2>0);
Sfilter2=Sfilter2-mean(Sfilter2(:));





%% recover alignments
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
lowResFluor2BF=alignments.lowResFluor2BF;
S2AHiRes=alignments.S2AHiRes;
Hi2LowResF=alignments.Hi2LowResF;
LowResBF2Hi=lowResFluor2BF;


temp1=invert(Hi2LowResF.t_concord);
temp2=(lowResFluor2BF.t_concord);
temp1=temp1.T*temp2.T;
LowResBF2Hi.t_concord.T=temp1;

rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;
%% calculate offset between frame position and zposition

zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%% set up low magvideos

aviFiles=dir([dataFolder filesep '20*.avi']);
aviFiles={aviFiles.name}';
aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
if length(aviFiles)==2
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behaviorMovie=[dataFolder filesep aviFiles{~aviFluorIdx}];
    fluorMovie=[dataFolder filesep aviFiles{aviFluorIdx}];
else
    display('Select avi files, behavior and then low mag fluor');
    movies=uipickfiles('FilterSpec',dataFolder);
    behaviorMovie=movies{1};
    fluorMovie=movies{2};
end

behaviorVidObj = VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);
bfSize=[behaviorVidObj.Height,behaviorVidObj.Width];
Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);



%% set up timing alignments and lookups
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');
bf2fluorLookup=interp1(fluorAll.frameTime,fluorIdxList,bfAll.frameTime,'linear');
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
hiResStackRange=1:max(hiResData.stackIdx);
hiResStackRange=hiResStackRange(~isnan(stack2BFidx));
BF2stackIdx=interp1(stack2BFidx(~isnan(stack2BFidx)),hiResStackRange,bfIdxList,'nearest');


%%
bins=8;
newZ2=bfAll.flashTrack;
newZ2(bfAll.flashLoc)=nan;
%newZ2=colNanFill(newZ2');
%newZ2=smooth(newZ2,5);
newZ2=(newZ2-nanmean(newZ2))/nanstd(newZ2);
newZ2=round(newZ2*2+bins/2);
newZ2(newZ2<0)=0;
newZ2(newZ2>bins)=bins;
newZ2=newZ2+1;
zList=unique(newZ2);
zList=zList(~isnan(zList));
%% calculate multiple backgrounds for BF
%progressbar(0);
medBfAll=nan(bfSize(1),bfSize(2),length(zList));
meanBfAll=medBfAll;
for iZ=1:length(zList);
    timeList=find(newZ2==zList(iZ));
    fluor=zeros(bfSize(1),bfSize(2));
    
    for iTime=1:length(timeList);
        progressbar(iTime/length(timeList),iZ/max(newZ2));
        currentTime=timeList(iTime);
        bfFrame = read(behaviorVidObj,currentTime);
        fluor=fluor+double(bfFrame(:,:,1));
    end
    
    %  medBfAll(:,:,iZ)=median(bfStack,3);
    meanBfAll(:,:,iZ)=(fluor/length(timeList));
    
    
    %  progressbar(iZ/max(newZ2));
    
end

%% calculate  background for fluor
progressbar(0);
fluorStack=0;
skip=10;
counter=0;
for iTime=1:skip:length(fluorAll.frameTime);
    progressbar(iTime/length(fluorAll.frameTime));
    if ~any(fluorAll.flashLoc==iTime)
        fluorFrame = read(fluorVidObj,iTime);
        fluorStack=fluorStack+double(fluorFrame(:,:,1));
        counter=counter+1;
    end
end
%  medBfAll(:,:,iZ)=median(bfStack,3);
fluorFrameProj=(fluorStack/counter);
%%
imagesc(fluorFrameProj);
hold on
headRect=roipoly();
fluorBackground=fluorFrameProj;
fluorBackground(headRect)=nan;
fluorBackground=inpaint_nans(fluorBackground);
imagesc(fluorBackground);
%%

imagesc(mean(meanBfAll,3));
headRect=roipoly();
meanBfAll2=meanBfAll;

for iZ=1:size(meanBfAll,3);
    
    temp=meanBfAll2(:,:,iZ);
    temp(headRect)=nan;
    temp=inpaint_nans(temp);
    meanBfAll2(:,:,iZ)=temp;
end
imagesc(mean(meanBfAll2,3));
%% make filter background
meanBfAllFilt=meanBfAll;
for iZ=1:size(meanBfAll,3);
    temp=meanBfAll2(:,:,iZ);
    temp=bpass(temp,1,40);
    meanBfAllFilt(:,:,iZ)=temp;
end


%% View background subtraction
bfStack=[];
CLall=zeros(100,2,length(bfAll.frameTime));
%%
startFlag=1;

refIdx=12;
%%
close all
bfList=find(~isnan(bf2fluorLookup));
bfList=bfList(bfList<maxFrame);
nCells=16;
nFrames=length(bfList);
nSteps=ceil(nFrames/nCells);
bfCell=cell(nCells);
for i=1:nCells+1
    low=min((i-1)*nSteps+1,nFrames);
    high=min(nSteps*i,nFrames);
    bfFrameRaw = read(behaviorVidObj,bfList(low));
    display('Select Initial Points');
    fig=imagesc(bfFrameRaw);
    [xpts,ypts]=getpts();
    CL=[xpts,ypts];
    clStartI{i}=CL;
    pause(.3)
    if i<=nCells
        bfCell{i}=bfList(low:high);
        
    end
end
close all
bfCell=bfCell(1:end-1);
%%
bfCellRev=cellfun(@(x) flipud(x), bfCell, 'uniform',0);
bfCell=reshape([bfCell;bfCellRev],1,[]);
clStart=reshape([clStartI;circshift(clStartI,[0 -1])],1,[]);
%%
CLcell=cell(1,2*nCells);
CL_I=CLcell;
save([dataFolder filesep 'CLworkspace']);
%%
parfor iCell=1:nCells*2
    startFlag=1;
    iFrame=1;
    cellList=bfCell{iCell};
    
    CLall=zeros(100,2,length(cellList));
    IsAll=zeros(100,length(cellList));
    cm_fluor=[26 26];
    behaviorVidObj = VideoReader(behaviorMovie);
    fluorVidObj= VideoReader(fluorMovie);
    %%
    for iFrame=1:length(cellList);
        %%
        iTime=bfCell{iCell}(iFrame);
        fluorTime=round(bf2fluorLookup(iTime));
        tic
        try
            if ~isnan(newZ2(iTime)) && ~any(fluorAll.flashLoc==fluorTime) && fluorTime>=1
                %%
                
                bfFrameRaw = read(behaviorVidObj,iTime);
                bfFrameRaw=double(bfFrameRaw(:,:,1));
                backGround=meanBfAllFilt(:,:,newZ2(iTime));
                
                
                fluorFrameRaw=read(fluorVidObj,fluorTime);
                fluorFrameRaw=double(fluorFrameRaw(:,:,1))-fluorBackground;
                fluorFrameRaw(fluorFrameRaw<0)=0;
                fluorFrame2=fluorFrameRaw;
                % fluorFrame2=imwarp(fluorFrameRaw,lowResFluor2BF.t_concord,...
                %    'OutputView',lowResFluor2BF.Rsegment);
                
                bfFrameRaw2=bpass((bfFrameRaw),3,40);
                fluorFrame2=bpass(fluorFrame2,1,40);
                %additional scaling factor
                scalingCorr=sum(sum(backGround.*bfFrameRaw2))/sum(backGround(:).^2);
                bfFrame=bfFrameRaw2-backGround*scalingCorr*.9;
                bfFrame(bfFrame<0)=0;
                %calculate centroid of fluor image
                fluorMask=false(size(fluorFrame2));
                fluorMask(round(cm_fluor(2))+(-25:25),round(cm_fluor(1))+(-25:25))=true;
                fluorBW=(fluorFrame2>(max(fluorFrame2(:))/5));
                fluorFrameMask=fluorBW.*fluorMask;
                if (any(fluorFrameMask(:))) && ~startFlag
                    fluorBW=fluorFrameMask;
                end
                
                fluorBW=AreaFilter(fluorBW); %select only largest obj
                
                [cmy,cmx]=find(fluorBW);
                cm_fluor=mean([cmx cmy]);
                cm=transformPointsForward(lowResFluor2BF.t_concord,cm_fluor);
                inputImage=filter2(gaussFilter,bfFrame,'same');
                inputImage=inputImage-min(inputImage(:));
                inputImage=normalizeRange(inputImage);
                maxThresh=quantile(inputImage(inputImage>.2),.90);
                inputImage(inputImage>maxThresh)=maxThresh;
                inputImage=inputImage/maxThresh;
                %%
                if startFlag
                    
                    CLold=clStart{iCell};
                    CLold=distanceInterp(CLold,100);
                    
                    
                    [CLOut,Is]=ActiveContourFit_wormRef2(inputImage, cline_para, CLold,refIdx,cm);
                    startFlag=0;
                else
                    
                    oldTime=iFrame-1;
                    while isnan(newZ2(bfCell{iCell}(oldTime)));
                        oldTime=oldTime-1;
                    end
                    %                 CLHistory_x=squeeze(CLall(:,1,oldTime-2:oldTime));
                    %                 CLHistory_y=squeeze(CLall(:,2,oldTime-2:oldTime));
                    %                 CL_predict_x=CLHistory_x*[1 -3 3]';
                    %                 CL_predict_y=CLHistory_y*[1 -3 3]';
                    %               CLold=[CL_predict_x CL_predict_y];
                    CLold=CLall(:,:,oldTime);
                    [CLOut,Is,Eout]=ActiveContourFit_wormRef2(inputImage, cline_para, CLold,refIdx,cm);
                    %%
                    % if there's a kink at the head, do something about it;
                    kink=diff(sqrt(sum(bsxfun(@minus,CLOut([refIdx 2*refIdx],:),CLOut(1,:)).^2,2)))
                    if kink<-10 || doesCross(CLOut(1:30,:));
                        cline_para2=cline_para;
                        keyboard
                        cline_para2.gamma=cline_para.gamma;
                        cline_para2.kappa=20;
                        cline_para2.refSpring=.01;
                        CLold=CLold(refIdx:end,:);
                        CLold=distanceInterp(CLold,100);
                        
                        [CLOut,Is,Eout]=ActiveContourFit_wormRef2(inputImage,...
                            cline_para, CLold,refIdx,cm);
                        [bool, x, y]=doesCross(CLOut(1:30,:));
                        if bool
                            CLOut(x:(1+y),:)=[];
                        CLOut=distanceInterp(CLOut,100);
                        end
                        
                    end
                end
                %%
                %             s=sqrt(sum(diff(CLOut).^2,2));
                %             meanS=mean(s);
                %             if isnan(meanS)
                %                 meanS=6;
                %             end
                
                %    cline_para.refL=median([6,meanS,6]);
                if ~mod(iTime,show2)
                    imagesc(inputImage);
                    % colormap gray
                    hold on
                    plot(CLOut(:,1),CLOut(:,2),'r');
                    plot(CLOut([1 end],1),CLOut([1 end],2),'og');
                    scatter(cm(1),cm(2),'gx');
                    plot([CLOut(refIdx,1) cm(1)],[CLOut(refIdx,2) cm(2)],'g');
                    %quiver(xyzs(:,1),xyzs(:,2),fs(:,1),fs(:,2),'black');
                    hold off
                    drawnow
                end
                CLOut=distanceInterp(CLOut,100);
                
                CLall(:,:,iFrame)=CLOut;
                IsAll(:,iFrame)=Is;
            else
                CLall(:,:,iFrame)=CLall(:,:,iFrame-1);
                IsAll(:,iFrame)=IsAll(:,iFrame-1);
            end
            
            
            
            display(['Completed frame ' num2str(iTime) ' in ' num2str(toc) ' s'])
        catch me
            display(['error frame ' num2str(iTime) ', cell ' num2str(iCell)])
            rethrow(me)
        end
    end
    CLcell{iCell}=CLall;
    CL_I{iCell}=IsAll;
end



%% reorganize outputs
for iCell=2:2:nCells*2
    CLcell{iCell}=flip(CLcell{iCell},3);

    CL_I{iCell}=flip(CL_I{iCell},2);
end



%%
CL_Lenght_Fun=@(x)  squeeze(sum(sqrt(sum(diff(x).^2,2))))';

CLcell=reshape(CLcell,2,nCells);
CL_I=reshape(CL_I,2,nCells);
CL_Iall=cell2mat(CL_I(:)');
%%
CL_bright_check=mean(CL_Iall>graythresh(CL_Iall));
CL_Iavg=trimmean(CL_Iall(:,CL_bright_check>.4),30,2);

CL_Isum=cellfun(@(x) sum(bsxfun(@minus,x,CL_Iavg).^2), CL_I,'uniform',0);
CL_Csum=cellfun(@(x) squeeze(sum(sqrt(sum((diff(x,2)).^2,2)),1))',CLcell,'uniform',0);
CL_Isum=cell2mat(CL_Isum);
CL_Csum=cell2mat(CL_Csum);
CL_Isum=CL_Isum+CL_Csum/2;
CL_Isum(1,:)=smooth(CL_Isum(1,:),300);
CL_Isum(2,:)=smooth(CL_Isum(2,:),300);
[~,idx]=min(CL_Isum);

CL_length=cellfun(@(x) CL_Lenght_Fun(x),CLcell,'uniform',0);
CL_length=cell2mat(CL_length);
% centerline lengths for forward fitting and backward fitting
CL_f=cell2mat(permute(CLcell(1,:),[1 3 2]));
CL_b=cell2mat(permute(CLcell(2,:),[1 3 2]));
CL_If=cell2mat(CL_I(1,:));
CL_Ib=cell2mat(CL_I(2,:));
CL_Iout=CL_If;
CL=CL_f;
CL(:,:,idx==2)=CL_b(:,:,idx==2);
CL_Iout(:,idx==2)=CL_Ib(:,idx==2);

%%
centerline=zeros(size(CL,1),2,length(bfAll.frameTime));
centerline(:,:,bfList)=CL;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jn3'],'centerline');
centerline=zeros(size(CL,1),2,length(bfAll.frameTime));
centerline(:,:,bfList)=CL_f;
centerline=flip(centerline,2);

save([dataFolder filesep 'centerline_jnf'],'centerline');
centerline=zeros(size(CL,1),2,length(bfAll.frameTime));
centerline(:,:,bfList)=CL_b;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jnb'],'centerline');

