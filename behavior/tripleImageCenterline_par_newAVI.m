
maxFrame=Inf; %%% CHECK FOR LAST USEFUL FRAME IN VIDEO

cline_para.refIdx=10;

cline_para.tipRegion=45;
cline_para.endRepulsion=.3;
cline_para.repulsionD=20;
cline_para.heat=3;
cline_para.CLalpha=5;
cline_para.CLbeta=100;
cline_para.gamma=25;
cline_para.kappa=30;
cline_para.endkappa=20;
cline_para.gradient_force=60;
cline_para.showFlag=0;
cline_para.iterations=400;

cline_para.stretching_force_factor=[.3 .3];
cline_para.refSpring=.01;
cline_para.stretch_ends_flag=1;
cline_para.refL=6;
cline_para.memForce=.005;
close all
gaussFilter=fspecial('gaussian',30,5);%fspecial('gaussian',10,75);
gaussFilter2=fspecial('gaussian',50,15);%fspecial('gaussian',10,75);
%%
dataFolder=uipickfiles();
dataFolder=dataFolder{1};


%% set up different kernals

gaussKernal2=gausswin(200);
gaussKernal2=convnfft(gaussKernal2,gaussKernal2');



Sfilter=max(gaussKernal2(:))-gaussKernal2;



%% load alignment data
try
    %try to load the alignment file in the folder, otherwise, select them
    %individual in the registration folder
    alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
catch
display('Select Low Res Alignment')

lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);

display('Select Hi to Low Fluor Res Alignment')
Hi2LowResF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowResF=load(Hi2LowResF{1});

 display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;


%% if you select them individually, bundle them and save it for later use
alignments.lowResFluor2BF=lowResFluor2BF;
alignments.S2AHiRes=S2AHiRes;
alignments.Hi2LowResF=Hi2LowResF;

save([dataFolder filesep 'alignments'],'alignments');

end


%% set up low magvideos

camData=importdata([ dataFolder  filesep 'camData.txt']);
time=camData.data(:,2);
fluorMovie=[dataFolder filesep 'cam0.avi'];
behaviorMovie=[dataFolder filesep 'cam1.avi'];
NFrames=length(camData.data);


behaviorVidObj = VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);
bfSize=[behaviorVidObj.Height,behaviorVidObj.Width];
refPointsx=meshgrid(1:20:1088);
refPointsy=refPointsx';
refPoints=sub2ind([1088,1088],refPointsx(:),refPointsy(:));
refIntensity=nan(length(refPoints),NFrames);
parfor_progress(NFrames);
parfor iTime=1:NFrames;
   % if ~any(fluorAll.flashLoc==iTime)
   
   behaviorVidObj = VideoReader(behaviorMovie);
        bfFrame =read(behaviorVidObj,[iTime]);
        refIntensity(:,iTime)=bfFrame(refPoints);
        parfor_progress;
 %   end
end

%%
behaviorVidObj = VideoReader(behaviorMovie);
refIntensityM=mean(refIntensity);
refIntensityM=refIntensityM-mean(refIntensityM);
newZ2=refIntensityM;
flashLoc=newZ2>(4*std(newZ2));
refIntensityM=refIntensityM-mean(refIntensityM(:,~flashLoc));


refIntensityMean=mean(refIntensity(:,~flashLoc),2);
refIntensityZ=bsxfun(@minus, double(refIntensity),refIntensityMean);
refIntensityZ(:,flashLoc)=nan;
[coeff,score,latent,tsquared] = pca(refIntensityZ');
flashLoc=find(flashLoc);


%%
bins=8;
newZ2=colNanFill(score(:,1));
gaussKernal=gausswin(1000, 4);
gaussKernal=gaussKernal/sum(gaussKernal);
newZ2=newZ2-conv(newZ2,gaussKernal,'same');
newZ2(flashLoc)=nan;
%newZ2=colNanFill(newZ2');
%newZ2=smooth(newZ2,5);
boundaryI=quantile(newZ2,10);
newZ2=sum(bsxfun(@le,newZ2,boundaryI),2);
newZ2=newZ2+1;
zList=unique(newZ2);
zList=zList(~isnan(zList));
%% calculate multiple backgrounds for BF
%progressbar(0);
medBfAll=nan(bfSize(1),bfSize(2),length(zList));
meanBfAll=medBfAll;
parfor_progress(NFrames);

parfor iZ=1:length(zList);
    timeList=find(newZ2==zList(iZ));
    fluor=zeros(bfSize(1),bfSize(2));
    behaviorVidObj = VideoReader(behaviorMovie);

    for iTime=1:length(timeList);
        currentTime=timeList(iTime);
        bfFrame = read(behaviorVidObj,currentTime);
        fluor=fluor+double(bfFrame(:,:,1));
        parfor_progress
    end
    
    %  medBfAll(:,:,iZ)=median(bfStack,3);
    meanBfAll(:,:,iZ)=(fluor/length(timeList));
    
    
    %  progressbar(iZ/max(newZ2));
    
end
    behaviorVidObj = VideoReader(behaviorMovie);

%% calculate  background for fluor
progressbar(0);
fluorStack=0;
skip=10;
counter=0;
for iTime=1:skip:NFrames;
    progressbar(iTime/NFrames);
    if ~any(flashLoc==iTime)
        fluorFrame = read(fluorVidObj,iTime);
        fluorStack=fluorStack+double(fluorFrame(:,:,1));
        counter=counter+1;
    end
end
progressbar(1);

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
meanBfBW=meanBfAll2>90; %flip signs in bright area 
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
CLall=zeros(100,2,NFrames);
%%
startFlag=1;

%%
close all

bfList=1:NFrames;
nCells=16;
nFrames=length(bfList);
nSteps=ceil(nFrames/nCells);
bfCell_i=cell(nCells,1);
for i=1:nCells+1
    low=min((i-1)*nSteps+1,nFrames);
    high=min(nSteps*i,nFrames);
    bfFrameRaw = read(behaviorVidObj,bfList(low));
    display(['Select Points for frame ' num2str(i) ' of ' num2str(nCells+1)]);
    fig=imagesc(bfFrameRaw);
    [xpts,ypts]=getpts();
    CL=[xpts,ypts];
    clStartI{i}=CL;
    pause(.3)
    if i<=nCells
        bfCell_i{i}=bfList(low:high);
        
    end
end
close all

 %% preview centerlines
 close all
for i=1:nCells+1
    low=min((i-1)*nSteps+1,nFrames);
    bfFrameRaw = read(behaviorVidObj,bfList(low));
    imagesc(bfFrameRaw)
    hold on
    plot(clStartI{i}(:,1),clStartI{i}(:,2));
    pause(.3)
hold off
end

%bfCell_i=bfCell_i(1:end-1);
%%
bfCellRev=cellfun(@(x) fliplr(x), bfCell_i, 'uniform',0);
bfCell=reshape([bfCell_i,bfCellRev]',1,[]);
clStart=reshape([clStartI;circshift(clStartI,[0 -1])],1,[]);
%%
CLcell=cell(1,2*nCells);
CL_I=CLcell;
save([dataFolder filesep 'CLworkspace']);
%%
parfor_progress(nFrames);
show2=0;
tstart=tic;
%%
for iCell=1:2*nCells
    startFlag=1;
    iFrame=1;
    cellList=bfCell{iCell};
    CLall=CLcell{iCell};
    if isempty(CLall)
        CLall=zeros(100,2,length(cellList));
    end
    IsAll=CL_I{iCell};
    if isempty(IsAll)
        IsAll=zeros(100,length(cellList));
    end

    cm_fluor=[26 26];
    behaviorVidObj = VideoReader(behaviorMovie);
    fluorVidObj= VideoReader(fluorMovie);
    %%
    for iFrame=i:length(cellList);
        %%
        iTime=bfCell{iCell}(iFrame);
        tic
        try
            if ~isnan(newZ2(iTime)) && ~any(flashLoc==iTime) && iTime>=1
                    %%
                
                bfFrameRaw = read(behaviorVidObj,iTime,'native');
                bfFrameRaw=double(bfFrameRaw.cdata);
                backGroundRaw=meanBfAll2(:,:,newZ2(iTime));
                c=sum(sum(bfFrameRaw.*backGroundRaw))/sum(sum(backGroundRaw.^2));
                bfFrameRaw2=bfFrameRaw-backGroundRaw*c;
                %bfFrameRaw2=abs(bfFrameRaw2);;
                bfFrame=bfFrameRaw2;
                bfFrame=normalizeRange(bpass(bfFrame,2,80));
                bfFramestd=stdfilt(bfFrameRaw2,getnhood(strel('disk',5)));
                bfFramestd=bpass((bfFramestd),1,80);

                fluorFrameRaw=read(fluorVidObj,iTime,'native');
                fluorFrameRaw=double(fluorFrameRaw.cdata)-fluorBackground;
                fluorFrameRaw(fluorFrameRaw<0)=0;
                fluorFrame2=fluorFrameRaw;
                fluorFrame2=bpass(fluorFrame2,1,40);
                
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
                inputImage=2*bfFrame/max(bfFrame(:))+normalizeRange(bfFramestd);
                inputImage(inputImage<0)=0;
                %       inputImage=sqrt(inputImage);
                inputImage=filter2(gaussFilter,inputImage,'same');
                inputImage=inputImage-min(inputImage(:));
                inputImage=normalizeRange(inputImage);
                maxThresh=quantile(inputImage(inputImage>.2),.90);
                inputImage=inputImage/maxThresh;
                bfFramestd=normalizeRange(imfilter(bfFramestd,gaussFilter2));
                nnz(bfFrame>.1)
                %%
                if startFlag
                    
                    CLold=clStart{iCell};
                    CLold=distanceInterp(CLold,100);
                    
                    
                    [CLOut,Is]=ActiveContourFit_wormRef4(inputImage,bfFramestd, cline_para, CLold,refIdx,cm);
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
                    [CLOut,Is,Eout]=ActiveContourFit_wormRef4(inputImage,bfFramestd, cline_para, CLold,refIdx,cm);
                    %%
                    % if there's a kink at the head, do something about it;
                    kink=diff(sqrt(sum(bsxfun(@minus,CLOut([refIdx 2*refIdx],:),CLOut(1,:)).^2,2)));
                    if kink<-10 || doesCross(CLOut(1:30,:));
                        cline_para2=cline_para;
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
                    subplot(1,2,1)
                    imagesc(bfFrameRaw2);
                    % colormap gray
                    hold on
                    plot(CLOut(:,1),CLOut(:,2),'r');
                    plot(CLOut([1 end],1),CLOut([1 end],2),'og');
                    scatter(cm(1),cm(2),'gx');
                    plot([CLOut(refIdx,1) cm(1)],[CLOut(refIdx,2) cm(2)],'g');
                    %quiver(xyzs(:,1),xyzs(:,2),fs(:,1),fs(:,2),'black');
                    hold off
                    subplot(1,2,2);
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
            
            
                  %  parfor_progress;

            display(['Completed frame '  num2str(iTime) ', cell ' num2str(iCell) ...
                ' in ' num2str(toc) ' s'])
        catch me
                 %   parfor_progress;

            display(['error frame ' num2str(iTime) ', cell ' num2str(iCell)])
            
            rethrow(me)
        end
    end
    CLcell{iCell}=CLall;
    CL_I{iCell}=IsAll;
end
toc(tstart)



%%


%%
CL_Lenght_Fun=@(x)  squeeze(sum(sqrt(sum(diff(x).^2,2))))';

CLcell_2=reshape(CLcell,2,nCells);
CL_I_2=reshape(CL_I,2,nCells);

CLcell_2(2,:)=cellfun(@(x) flipdim(x,3),CLcell_2(2,:),'uniform',0);
CL_I_2(2,:)=cellfun(@(x) flipdim(x,3),CL_I_2(2,:),'uniform',0);

CL_Iall=cell2mat(CL_I_2(:)');

%%
CL_bright_check=mean(CL_Iall>graythresh(CL_Iall));
CL_Iavg=trimmean(CL_Iall(:,CL_bright_check>.4),30,2);

CL_Isum=cellfun(@(x) sum(bsxfun(@minus,x,CL_Iavg).^2), CL_I_2,'uniform',0);
CL_Csum=cellfun(@(x) squeeze(sum(sqrt(sum((diff(x,2)).^2,2)),1))',CLcell_2,'uniform',0);
CL_Isum=cell2mat(CL_Isum);
CL_Csum=cell2mat(CL_Csum);
CL_Isum=CL_Isum+CL_Csum/2;
CL_Isum(1,:)=smooth(CL_Isum(1,:),300);
CL_Isum(2,:)=smooth(CL_Isum(2,:),300);
[~,idx]=min(CL_Isum);

CL_length=cellfun(@(x) CL_Lenght_Fun(x),CLcell_2,'uniform',0);
CL_length=cell2mat(CL_length);
% centerline lengths for forward fitting and backward fitting
CL_f=cell2mat(permute(CLcell_2(1,:),[1 3 2]));
CL_b=cell2mat(permute(CLcell_2(2,:),[1 3 2]));
CL_If=cell2mat(CL_I_2(1,:));
CL_Ib=cell2mat(CL_I_2(2,:));
CL_Iout=CL_If;
CL=CL_f;
CL(:,:,idx==2)=CL_b(:,:,idx==2);
CL_Iout(:,idx==2)=CL_Ib(:,idx==2);

%%
centerline=zeros(size(CL,1),2,NFrames);
centerline(:,:,bfList)=CL;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jn3'],'centerline');
centerline=zeros(size(CL,1),2,NFrames);
centerline(:,:,bfList)=CL_f;
centerline=flip(centerline,2);

save([dataFolder filesep 'centerline_jnf'],'centerline');
centerline=zeros(size(CL,1),2,NFrames);
centerline(:,:,bfList)=CL_b;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jnb'],'centerline');

