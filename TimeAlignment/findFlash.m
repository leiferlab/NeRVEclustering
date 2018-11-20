function  imFlash=findFlash(imFolderIn,flag)
% function takes folder of tif images and searches for the flashes in the
% images by finding the average intensities in either the entire image or
% some region of the image.
mostRecent=getappdata(0,'mostRecent');
if nargin==0
    imFolderIn=uipickfiles('FilterSpec',mostRecent);
    setappdata(0,'mostRecent',imFolderIn{1});
    %  imFolder=imFolder{1};
    flag.customRoi=0;
    flag.sparsesearch=0;
    flag.parPool=0;
end

if nargin==1
    flag.customRoi=0;
    flag.sparsesearch=0;
    flag.parPool=0;
end
if ~iscell(imFolderIn)
    imFolderIn={imFolderIn};
end
if ~isfield(flag, 'custormRoi')
    flag.customRoi=0;
end
if ~isfield(flag, 'sparsesearch')
    flag.sparseserach=0;
end

%% load images and find flash in images using user defined ROI
for iFolder=1:length(imFolderIn)
    imFolder=imFolderIn{iFolder};
    if isdir(imFolder)
        imFiles=dir([imFolder filesep '*.tif']);
        vidObj=[];
        stackSize=length(imFiles);
        initialIm=(imread([imFolder filesep imFiles(1).name], 'tif'));
        aviFlag=0;
    elseif strfind(imFolder, '.avi')
        vidObj = VideoReader(imFolder);
        initialIm = sum(read(vidObj, inf),3);
        stackSize= vidObj.NumberOfFrames;
        aviFlag=1;
        imFiles=zeros(1,stackSize);
    end
    
    if flag.customRoi
        progressbar;
        for i=round(stackSize/2):round(stackSize/2)+200;
            if aviFlag;
                temp= read(vidObj, i);
                temp=(sum(double(temp),3));
            else
                temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
            end
            initialIm=max(initialIm,temp);
            progressbar((i-round(stackSize/2))/200);
            
        end
        imsize=size(initialIm);
        fig=imagesc(initialIm);
        display('Select area to find flash');
        roiFlash=roipoly;
        delete(fig)
    else
        roiFlash=true(size(initialIm));
    end
    
    % fig=imagesc(initialIm);
    % display('Select an cell to try to align');
    % roiCell1=roipoly;
    % delete(fig)
    %
    % fig=imagesc(initialIm);
    % display('Select another');
    % roiCell2=roipoly;
    % close all
    
    
    %% search for the flash
    imFlash=zeros(1,stackSize);
    %cell1=imFlash;
    %cell2=imFlash;
    if flag.sparsesearch
        progressbar;
        for i=1:10:stackSize/10
            if aviFlag;
                temp= read(vidObj, i);
                temp=(sum(double(temp),3));
            else
                temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
            end
            imFlash(i)=mean(temp(roiFlash));
            %     cell1(i)=mean(temp(roiCell1));
            %     cell2(i)=mean(temp(roiCell2));
            progressbar(i/stackSize*10);
        end
        progressbar(1);
        flashMean=mean(imFlash(imFlash>0));
        flashStd=std(imFlash(imFlash>0));
        
        
        while ~any(imFlash>(flashMean+5*flashStd))
            i=i+10;
            if aviFlag;
                temp= read(vidObj, i);
                temp=(sum(double(temp),3));
            else
                temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
            end
            imFlash(i)=mean(temp(roiFlash));
            flashMean=mean(imFlash(imFlash>0));
            flashStd=std(imFlash(imFlash>0));
            
        end
        i=find(imFlash==max(imFlash));
        for k=max(i-100,0):min(i+100,stackSize)
            if aviFlag;
                temp= read(vidObj, i);
                temp=(sum(double(temp),3));
            else
                temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
            end
            imFlash(k)=mean(temp(roiFlash));
            progressbar((k-i+100)/200);
            
            
        end
        imFlash(imFlash==0)=flashMean;
        
    elseif flag.parPool
        parfor_progress(stackSize)
        parfor i=1:stackSize
            if aviFlag;
                vidObjpar = VideoReader(imFolder);
                
                temp= read(vidObjpar, i);
                temp=temp(:,:,1);
            else
                temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
            end
            imFlash(i)=mean(temp(roiFlash));
            %     cell1(i)=mean(temp(roiCell1));
            %     cell2(i)=mean(temp(roiCell2));
            %progressbar(i/stackSize);
            parfor_progress;
        end
        parfor_progress(0);
    else
        for i=1:stackSize
            if aviFlag;
                temp= read(vidObj, i);
                temp=temp(:,:,1);
            else
                temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
            end
            imFlash(i)=mean(temp(roiFlash));
            %     cell1(i)=mean(temp(roiCell1));
            %     cell2(i)=mean(temp(roiCell2));
            %progressbar(i/stackSize);
            progressbar(i/stackSize)
        end
        progressbar(1);
        
    end
    
    if aviFlag
        
        save(strrep(imFolder, '.avi','flashTrack'),'imFlash');
    else
        save([imFolder filesep 'flashTrack'],'imFlash');
    end
end
