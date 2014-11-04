

%load folder and extract and align timing data.

alignFlag=0;
tStack=1;
overwrite=1;
startPoint=1;


%load registration file
[regFile,regFolder]=uigetfile('Y:\CommunalCode\3dbrain\registration\');
load([regFolder filesep regFile]);



%%

%% segment subimages and create masks
[folderList] = uipickfiles; %select folders to analyze
se=strel('disk',5);
%%
%looping through folders to analyze
for iFolder=1:length(folderList)
    
    imFolder=folderList{iFolder};
    imNames=dir([imFolder filesep '*.tif']);
    mkdir([imFolder filesep 'stackData']);
    
    if ~overwrite
        matFiles=dir([imFolder filesep 'stackData' filesep '*.mat']);
        matFileNames={matFiles.name}';
        mat2tifFileNames=cellfun(@(x) strrep(x,'.mat','.tif'),matFileNames,...
            'UniformOutput',0);
        [Lia,Locb]=ismember(mat2tifFileNames,{imNames.name}');
        startPoint=max(max(Locb),startPoint);
    end
    
    
    maskStack=[];
    wormStack=[];
    activityStack=[];
    imNameStack=[];
    
    movieLength=length(imNames);
    %movieLength=2932;
    %% analyze all images
    
    for i=1:10;
        saveChunk=(movieLength/10);
        saveChunkI=round((i-1)*saveChunk+1:i*saveChunk);
        outputData=struct('centroids',[]);
        outputData=repmat(outputData,length(saveChunkI),1);
        
        
        parfor iImage=saveChunkI;
            tic
            outputFile=[imFolder filesep 'stackData' filesep  strrep(imNames(iImage).name,'.tif','.mat')];
            if ~exist(outputFile,'file')
                try
                    % load stacks
                    
                    temp=double(imread([imFolder filesep imNames(iImage).name],'tif'));
                    %temp=pixelIntensityCorrection(temp);
                    
                    % cut out gcamp and tdtomato sections of image
                    temp_activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
                    worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
                    %         temp_activity=temp((rect2(1)+1):rect2(3),(1+rect2(2)):rect2(4));
                    %         worm=temp((rect1(1)+1):rect1(3),(1+rect1(2)):rect1(4));
                    
                    temp_activity=imwarp(temp_activity,t_concord,'OutputView',Rsegment);
                    temp_activity(padRegion)=median(temp_activity(~padRegion));
                    %    activity=bpass_jn((temp_activity),1,[20,20]);
                    activity=temp_activity;
                    imsize=size(worm);
                    %do segmentation
                    %wormMask=WormSegmentHessian2D(worm);
                    wormMask=WormSegmentHessian2D_whole(worm);
                    wormMask= bwmorph(wormMask,'clean');
                    
                    
                    
                    %build stack
                    %     maskStack=  cat(3,maskStack,bwlabel(wormMask));
                    % wormStack= cat(3,wormStack,worm);
                    % activityStack=cat(3,activityStack,activity);
                    % imNameStack=cat(1,imNameStack,{ strrep(imNames(iImage).name,'.tif','.mat')});
                    
                    % for multiple images where we take a single good slice
                    %if size(maskStack,3)==tStack;
                    
                    
                    %pick the part of the stack you want
                    % wormMask=(mean(maskStack>0,3)>.8);
                    % worm=wormStack(:,:,ceil(tStack/2));
                    % activity=activityStack(:,:,ceil(tStack/2));
                    
                    
                    
                    %make labels and dilate, do the same in reverse order so theres no bias
                    wormLabelMask1=bwlabeln(wormMask);
                    wormLabelMask2=(max(wormLabelMask1(:))+1-wormLabelMask1).*(wormLabelMask1>0);
                    wormLabelMask1=imdilate(wormLabelMask1,se);
                    wormLabelMask2=imdilate(wormLabelMask2,se);
                    wormLabelMask2=(max(wormLabelMask2(:))-wormLabelMask2+1).*(wormLabelMask2>0);
                    
                    %calculate fluor intensities for each region of labeled masks
                    Rstats=regionprops(wormLabelMask1,worm,'Area','PixelValues','Centroid');
                    Gstats=regionprops(wormLabelMask1,activity,'PixelValues');
                    
                    
                    Rstats2=regionprops(wormLabelMask2,worm,'PixelValues');
                    Gstats2=regionprops(wormLabelMask2,activity,'PixelValues');
                    
                    Rintensities=cellfun(@(x) trimmean(x,10), {Rstats.PixelValues}');
                    Rintensities2=cellfun(@(x) trimmean(x,10), {Rstats2.PixelValues}');
                    
                    Gintensities=cellfun(@(x) trimmean(x,10), {Gstats.PixelValues}');
                    Gintensities2=cellfun(@(x) trimmean(x,10), {Gstats2.PixelValues}');
                    
                    Rintensities=.5*(Rintensities+Rintensities2);
                    Gintensities=.5*(Gintensities+Gintensities2);
                    
                    Volume=[Rstats.Area]';
                    centroids=cell2mat({Rstats.Centroid}');
                    %save outputs in unique file
                    
                    outputData(iImage).centroids=centroids;
                    outputData(iImage).Rintensities=Rintensities;
                    outputData(iImage).Gintensities=Gintensities;
                    outputData(iImage).Volume=Volume;
                    outputData(iImage).wormMask=uint8(wormLabelMask1);
                    
                    
                    
                    
                    display(['Completed stack' num2str(iImage,'%04d') ' in ' num2str(toc) ' seconds']);
                catch
                    display(['Error Caught in stack' num2str(iImage,'%04d') ' in ' num2str(toc) ' seconds']);
                end
                
            else
                display(['Skipping ' outputFile]);
            end
            
            
            
        end
        for iImage=saveChunkI
            outputFile=[imFolder filesep 'stackData' filesep  strrep(imNames(iImage).name,'.tif','.mat')];
            parsavestruct(outputFile,outputData(iImage));
        end
        
    end
end
