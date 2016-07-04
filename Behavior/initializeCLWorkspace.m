
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
%% if there's a background image, load it as well into alignments.
display('select a background image for this size himag video');

backgroundImage=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\background');
if iscell{backgroundImage}
backgroundImage=load(backgroundImage{1});
backgroundImage=backgroundImage.backgroundImage;
else
    backgroundImage=0;
end


%% if you select them individually, bundle them and save it for later use
alignments.lowResFluor2BF=lowResFluor2BF;
alignments.S2AHiRes=S2AHiRes;
alignments.Hi2LowResF=Hi2LowResF;
alignments.background=backgroundImage;



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

%% select centerlines

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


%% calculate backgrounds and projections

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


%bfCell_i=bfCell_i(1:end-1);
%%
bfCellRev=cellfun(@(x) fliplr(x), bfCell_i, 'uniform',0);
bfCell=reshape([bfCell_i,bfCellRev]',1,[]);
clStart=reshape([clStartI;circshift(clStartI,[0 -1])],1,[]);
%%
CLcell=cell(1,2*nCells);
CL_I=CLcell;
save([dataFolder filesep 'CLworkspace']);