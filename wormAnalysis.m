

%load folder and extract and align timing data.
display('Select Registration mat File');
[regFile,regFolder]=uigetfile('Y:\CommunalCode\3dbrain\registration\');

%%
display('Select folder to analyze');

imFolder=uigetdir;

if ~exist(fullfile(imFolder,'StackInfo.mat'),'file')
[initialIm,stackInfo]=timeTraceAnalysis(imFolder);
save([imFolder filesep 'stackInfo'],'stackInfo','initialIm');
end

%% load registration file
load([imFolder filesep 'stackInfo'],'stackInfo','initialIm');

    load([regFolder filesep regFile]);


%% segment subimages and create masks
paraPool=gcp('nocreate');
if ~isempty(paraPool)
paraPool=parpool('local',8,'IdleTimeout', 600);
else
    delete(paraPool)
    paraPool=parpool('local',8,'IdleTimeout', 600);

end
    
%%

options.thresh1=.03; %initial Threshold
options.hthresh=-.0001; %threshold for trace of hessian.
options.minObjSize=500; 
options.maxObjSize=Inf;
options.watershedFilter=1;
options.filterSize=[50,50,50];
options.pad=9;
options.noise=12;
options.show=0;
options.maxSplit=1;
options.minSphericity=.55;
%%


mkdir([imFolder filesep 'stackData']);


for iStack=1:length(stackInfo);
    tic
stackSize=length(stackInfo(iStack).fileNames);
    worm=zeros(rectSize1(2),rectSize1(1),stackSize);
    activity=worm;
    % load stacks
    for slice=1:stackSize
        
        temp=double(imread([imFolder filesep stackInfo(iStack).fileNames(slice).name],'tif'));
        temp=pixelIntensityCorrection(temp);
        temp_activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
        worm(:,:,slice)=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
        temp_activity=imwarp(temp_activity,t_concord,'OutputView',Rsegment);
        temp_activity(padRegion)=median(temp_activity(~padRegion));
        activity(:,:,slice)=temp_activity;
    end
    
    imsize=size(worm);  
    
    %resize image, arbitrary for now
   worm=image_resize(worm,imsize(1),imsize(2),2*imsize(3));
      activity=image_resize(activity,imsize(1),imsize(2),2*imsize(3));
%do segmentation
    wormMask=WormSegmentHessian3d(worm,options);
    
    %look up intensities on both channels
    wormLabelMask=bwlabeln(wormMask);
wormcc=bwconncomp(wormMask);
stats=regionprops(wormcc,'Centroid','Area');
centroids=reshape([stats.Centroid],3,[])';
Rintensities=cellfun(@(x) mean(worm(x)),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) mean(activity(x)),[wormcc.PixelIdxList])';

    %interpolate Z properly and scale
centroids(:,3)=interp1(stackInfo(iStack).z,centroids(:,3)/2)*100; %arb scaling for now
Volume=[stats.Area]';
realTime=interp1(stackInfo(iStack).time,centroids(:,3)/2); %arb scaling for now

%save outputs in unique file
outputFile=[imFolder filesep 'stackData' filesep 'stack' num2str(iStack,'%04d') 'data'];
outputData(iStack).centroids=centroids;
outputData(iStack).Rintensities=Rintensities;
outputData(iStack).Gintensities=Gintensities;
outputData(iStack).Volume=Volume;
outputData(iStack).time=realTime;
outputData(iStack).wormMask=wormLabelMask;

parsavestruct(outputFile,outputData(iStack));


display(['Completed stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);
end
