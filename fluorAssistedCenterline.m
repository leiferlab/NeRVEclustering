
%function [V,pointStats,Vproj,side,xyOffset2,wormBW2]=trippleImageCenterline(dataFolder,alignments,iStack,side,show)
%%

cline_para.CLalpha=5;
cline_para.CLbeta=50;
cline_para.gamma=10;
cline_para.kappa=20;
cline_para.gradient_force=50;
cline_para.showFlag=0;
cline_para.iterations=400;
cline_para.stretching_force_factor=[2 2];
cline_para.refSpring=0.5;
cline_para.stretch_ends_flag=1;
cline_para.refSpring=5;
cline_para.springForce=1;
cline_para.refL=6;
cline_para.memForce=.02;
close all
gaussFilter=fspecial('gaussian',20,75);
%%
try
    %% time align videos (using this also gets the unused high res videos)
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
    BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest');
    
    
    %% calculate multiple backgrounds for BF, depending on frames avg intensity
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
    %progressbar(0);
    medBfAll=nan(bfSize(1),bfSize(2),length(zList));
    meanBfAll=medBfAll;
    for iZ=1:length(zList);
        timeList=find(newZ2==zList(iZ));
        fluor=nan(bfSize(1),bfSize(2),length(timeList));
        
        for iTime=1:length(timeList);
            %    progressbar(iZ/max(newZ2),iTime/length(timeList));
            currentTime=timeList(iTime);
            bfFrame = read(behaviorVidObj,currentTime);
            fluor=fluor+double(bfFrame(:,:,1));
        end
        
        %  medBfAll(:,:,iZ)=median(bfStack,3);
        meanBfAll(:,:,iZ)=(fluor/length(timeList));
        
        
        progressbar(iZ/max(newZ2));
        
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
    tstart=max(iTime-50,1);
    for iTime=tstart:1:150000
        fluorTime=round(bf2fluorLookup(iTime));
        if ~isnan(newZ2(iTime)) && ~isnan(iStack) && ~any(fluorAll.flashLoc==fluorTime)
            %%
            tic
            bfFrameRaw = read(behaviorVidObj,iTime);
            bfFrameRaw=double(bfFrameRaw(:,:,1));
            backGround=meanBfAllFilt(:,:,newZ2(iTime));
            
            
            if isnan(fluorTime) |fluorTime<1
                fluorTime=1;
            end
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
            inputImage=normalizeRange(inputImage);
            %%
            if startFlag
                display('Select Initial Points');
                fig=imagesc(bfFrame);
                [xpts,ypts]=getpts();
                CL=[xpts,ypts];
                delete(fig);
                CL=distanceInterp(CL,100);
                
                
                CLOut=ActiveContourFit_wormRef(inputImage, cline_para, CL,refIdx,cm);
                startFlag=0;
            else
                
                oldTime=iTime-1;
                while isnan(newZ2(oldTime))
                    oldTime=oldTime-1;
                end
                %                 CLHistory_x=squeeze(CLall(:,1,oldTime-2:oldTime));
                %                 CLHistory_y=squeeze(CLall(:,2,oldTime-2:oldTime));
                %                 CL_predict_x=CLHistory_x*[1 -3 3]';
                %                 CL_predict_y=CLHistory_y*[1 -3 3]';
                %               CLold=[CL_predict_x CL_predict_y];
                CLold=CLall(:,:,oldTime);
                CLOut=ActiveContourFit_wormRef(inputImage, cline_para, CLold,refIdx,cm);
            end
            s=sqrt(sum(diff(CLOut).^2,2));
            meanS=mean(s);
            if isnan(meanS)
                meanS=6;
            end
            
            %    cline_para.refL=median([6,meanS,6]);
            if ~mod(iTime,1)
                imagesc(bfFrameRaw);
                % colormap gray
                hold on
                plot(CLOut(:,1),CLOut(:,2),'r');
                scatter(cm(1),cm(2),'gx');
                plot([CLOut(refIdx,1) cm(1)],[CLOut(refIdx,2) cm(2)],'g');
                %quiver(xyzs(:,1),xyzs(:,2),fs(:,1),fs(:,2),'black');
                hold off
                drawnow
            end
            CLOut=distanceInterp(CLOut,100);
            
        end
        CLall(:,:,iTime)=CLOut;
        
        
        display(['Completed frame ' num2str(iTime) ' in ' num2str(toc) ' s'])
    end
    
    
    %%
    %%
    
    fclose(Fid);
catch me
    fileName=[imageFolder2 filesep 'ERROR' num2str(iStack,'%3.5d')];
    save(fileName);
    rethrow(me)
end
