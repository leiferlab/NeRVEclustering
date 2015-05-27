imFolder='F:\20141212\BrainScanner20141212_145951\CLstraight2\';

zindexer=@(x,s) x./(s)+1;

options.thresh1=0.05;
options.minObjSize=50;
options.filterSize=[10,10,4];
imageRange=500:700;
%%
for i=1:length(imageRange)
    %%
    iImage=imageRange(i);
    imFile=[imFolder filesep 'image00' num2str(iImage) '.tif'];
    mapFile=[imFolder filesep 'imageMap00' num2str(iImage) '.tif'];
    controlFile=[imFolder filesep 'controlPoints' num2str(iImage,'%3.5d')];
    
    controlPoints=load(controlFile,'Fpoints2');
    controlPoints=controlPoints.Fpoints2;
    image=stackLoad(imFile);
    image(isnan(image))=0;
    
    mapImage=stackLoad(mapFile);
    if mod(iImage,2)==0;
        image=flip(image,3);
        mapImage=flip(mapImage,3);
        controlPoints(:,3)=201-controlPoints(:,3);
        
    end
    
    imsize=size(image);
    [wormBW2,wormtop]=WormSegmentHessian3dStraighten(image,options);
    
    BWplot=(squeeze(sum(sum(imdilate(wormBW2,true(10,10,10)),1),2)));
    BWplot=smooth(BWplot,20);
    [~,locs]=findpeaks(BWplot);
    endpts=locs([1,end]);
    [~,locs]=findpeaks(-BWplot);
    botpoint1=locs((locs>endpts(1)));
    if isempty(botpoint1);botpoint1=1;end;
    botpoint1=botpoint1(1);
    botpoint2=locs((locs<endpts(2)));
    if isempty(botpoint2);botpoint2=imsize(3);end;
    
    botpoint2=botpoint2(end);
    botpoint1(botpoint1>imsize(3)*1/4)=1;
    botpoint2(botpoint2<imsize(3)*3/4)=imsize(3);
    
    
    cc=bwconncomp(wormBW2);
    
    badRegions=(cellfun(@(x) any(zindexer(x,imsize(1)*imsize(2))<=botpoint1),cc.PixelIdxList)...
        |cellfun(@(x) any(zindexer(x,imsize(1)*imsize(2))>=botpoint2),cc.PixelIdxList))';
    
    wormBW2(cell2mat(cc.PixelIdxList(badRegions)'))=false;
    cc.PixelIdxList=cc.PixelIdxList(~badRegions);
    cc.NumObjects=nnz(~badRegions);
    
    stats=regionprops(cc,image,'Centroid','MeanIntensity',...
        'Area');
    
    statsMap=regionprops(cc,mapImage,'MeanIntensity');
    regionLabel=round([statsMap.MeanIntensity]);
    intensities=[stats.MeanIntensity]';
    P=[cell2mat({stats.Centroid}'),regionLabel',iImage*ones(cc.NumObjects,1)...
        (1:cc.NumObjects)'  intensities];
    P(:,[1 2])=P(:,[2 1]);
    options.method='invdist';
    options.radius=20;
    options.power=1;
    
    
    minLength=min(length(controlPoints),length(refPoints));
    noNans=~any(isnan([controlPoints(1:minLength,:), refPoints(1:minLength,:)]),2);
    moving=controlPoints(noNans,:);
    model=refPoints(noNans,:);
    fiducialsAll(i).moving=moving;
    fiducialsAll(i).model=model;
    
    %%
    if i>1
        [imOut,Xw,Yw,Zw]  = tpswarp3(image, size(image),...
            moving(:,[2 1 3]),model(:,[2 1 3]));
        %      controlmoved=tpswarp3points(controlAll{i}(1:minLength,:),refPoints(1:minLength,:),controlAll{i}(1:minLength,:));
        controlmoved=tpswarp3points(moving(:,[2 1 3]),model(:,[2 1 3]),P(:,1:3));
    else
        controlmoved=P(:,1:3);
        refPoints=controlAll{1};
        imOut=image;
        
        refPointsLin=sub2ind(size(image),refPoints(:,2),refPoints(:,1),refPoints(:,3));
        refPosIdx=refPointsLin;
        refPosIdx(~isnan(refPointsLin))=mapImage(refPointsLin(~isnan(refPointsLin)));
        
    end
    fileName4=[imFolder filesep 'fiducialStraight' num2str(iStack,'%3.5d')];
    
    tiffwrite(fileName4,imOut,'tif');
    fiducialsAll(i).ID=refPosIdx(noNans);
    
    
    TrackData{i}=P;
    controlAllMoved{i}=controlmoved;
    display(['Completed ' num2str(iImage)]);
end
controlmoved=[];
model=[];
%%
param.dim=3;
param.excessive=4;
progressbar(0,0);
windowSearch=5;
for i=1:30%length(TrackData)
    outRange=1:30;%max(1,i-windowSearch):min(length(TrackData),i+windowSearch);
    TrackMatrixi=zeros(size(TrackData{i},1),length(outRange));
    
    for j=1:30%outRange;
        progressbar(i/length(TrackData),j/length(TrackData));
        [i j]
        for regionId=0:2
            
            
            
            T1=[TrackData{i} (1:length(TrackData{i}))'];
            T2=[TrackData{j} (1:length(TrackData{j}))'];
            select1=T1(:,4)==regionId & T1(:,end-1)>40;
            select2=T2(:,4)==regionId & T2(:,end-1)>40;
            T1=T1(select1,1:3);
            T2=T2(select2,1:3);
            if ~isempty(T2) && ~isempty(T1)
                if size(T1,1)>1 && size(T2,1)>1
                    [Transformed_M, multilevel_ctrl_pts, multilevel_param] = ...
                        gmmreg_L2_multilevel_jn(T2,T1, 3, [20, 5, 1], ...
                        [0.0008, 0.0000008, 0.00000008],[0 0],...
                        [0.00001 0.0001 0.001],0);
                else
                    Transformed_M=T2;
                    
                    
                end
                
                trackInput=[T1  T1 find(select1) ones(size(T1(:,1))); ...
                    Transformed_M T2 find(select2) 2*ones(size(Transformed_M(:,1)))];
                TrackOut=nan;
                counter=10;
                while(all(isnan(TrackOut(:))))
                    TrackOut=trackJN(trackInput,counter,param);
                    counter=counter-1;
                end
                
                TrackOut(:,1:3)=[];
                TrackStats=round(TrackOut(:,4:end));
                TrackedIDs=TrackStats([1;diff(TrackStats(:,3))]==0,end);
                TrackStats=TrackStats(ismember(TrackStats(:,end),TrackedIDs),:);
                track1=TrackStats(1:2:end,1);
                track2=TrackStats(2:2:end,1);
                TrackMatrixi(track1,j-outRange(1)+1)=track2;
            end
        end
    end
    TrackMatrix{i}=TrackMatrixi;
end

%% do matching of tracks by clustering similarities in matching matrices

indexAdd=cell2mat(cellfun(@(x) size(x,1),TrackMatrix,'uniform',0)');
indexAdd=[0; cumsum(indexAdd(1:end))]';

for i=1%:length(TrackMatrix)-windowSearch
    %%
    outRange=max(1,i-windowSearch):min(length(TrackData),i+windowSearch);
    for j=outRange
        outRange2{j}=max(1,j-windowSearch):min(length(TrackData),j+windowSearch);
    end
    %transitionMatrix=false(size(cell2mat(TrackMatrix(outRange2{1})),1));
    subTrackMatrix=TrackMatrix(outRange2{i});
    transitionMatrixSize=cellfun(@(x) size(x,1),subTrackMatrix,'uniform',0);
    transitionMatrixSize=sum(cell2mat(transitionMatrixSize));
    
    transitionMatrix=false(transitionMatrixSize);
    
    
    for j=outRange2{i}
        TrackMatrixi=TrackMatrix{j};
        [outRangeoverlap,ia, ib]=intersect(outRange2{j},outRange2{i});
        TrackMatrixi=TrackMatrixi(:,ia);
        TrackMatrixi(TrackMatrixi==0)=nan;
        TrackMatrixi(:,outRangeoverlap==j)=nan;
        startIdx=indexAdd(outRange(1));
        TrackMatrixi=bsxfun(@plus, TrackMatrixi, indexAdd(outRangeoverlap)-startIdx);
        startPos=indexAdd(j)-startIdx+(1:size(TrackMatrixi,1));
        startPos=repmat(startPos,length(outRangeoverlap),1)';
        
        
        transitionIdx=sub2ind(size(transitionMatrix),startPos(startPos(:)>0),TrackMatrixi(startPos(:)>0));
        transitionIdx=transitionIdx(~isnan(transitionIdx));
        transitionMatrix(transitionIdx)=true;
        
    end
    transitionMatrix=transitionMatrix+transitionMatrix'+eye(size(transitionMatrix));
    transitionMatrix=bsxfun(@rdivide, transitionMatrix,sum(transitionMatrix));
    transitionMatrix(isnan(transitionMatrix))=0;
    %%
    tcorr2=corr(transitionMatrix);
    
    
    corrmat=tcorr2.*~eye(size(tcorr2));
    corrmat=squareform(corrmat);
    Z=linkage(1-corrmat,'complete');
    c=cluster(Z,'cutoff',.5,'criterion','distance');
    
    caccum=accumarray(c,ones(size(c)));
    caccum=find(caccum);
    c(ismember(c,caccum))=0;
    [~,ia,ib]=unique(c);
    [~,ic]=sort(ia);
    [~, id]=sort(ic);
    c2=id(ib);
    c2(c==0)=0;
    ccell{i}=mat2cell(c2,diff(indexAdd(min(outRange):max(outRange+1))));
    
end

%%
ccell2=ccell;

%loop over all clusters
for i=1:2%:length(ccell2)-1
    outRange=outRange2{i};
    % outRange(outRange<=i)=[];
    %loop over window of cluster comprisons
    for j=outRange;
        assign1=ccell2{i};
        assign2=ccell2{j};
        

        [outRangeoverlap,ia, ib]=intersect(outRange2{i},outRange2{j});
        max1=max(unique(cell2mat(assign1)));
        assign1lin=cell2mat(assign1(ia));
        assign2lin=cell2mat(assign2(ib));
        reassign1= cell2mat(ccell2{i});

        reassign2= cell2mat(ccell2{j});
        reassign3=reassign2;
        maxVal=max(assign1lin);
        %loop over indeces
        for iIdx=1:max(assign2lin)
            newVal=assign1lin(assign2lin==iIdx);
            newVal=newVal(newVal~=0);
            [iIdx; newVal]
            if ~isempty(newVal)
                %    [i,iIdx, mode( newVal(newVal~=0))]
                 newVal(newVal~=0)
                
                newValMode=mode(newVal(newVal~=0));
                if newValMode~=0 && ~isempty(newValMode)
                    newVal=mode(newVal(newVal~=0));
                else
                    newVal=maxVal+1;
                    maxVal=newVal;
                end
            else
                newVal=maxVal+1;
                maxVal=newVal
            end
            
            newValAll(iIdx)=newVal;
            reassign3(reassign2==iIdx)=newVal;
            
        end
        reassign2=mat2cell(reassign3,diff(indexAdd(outRange2{j}(1)...
            :max(outRange2{j}(end)+1))));
        reassign2j=reassign2{outRange2{j}==i};
        reassign2i=assign1{outRange2{i}==i};
        reassign2i(reassign2i==0)=reassign2j(reassign2i==0);
        ccell2{i}{outRange2{i}==i}=reassign2i;
        ccell2{j}(outRange2{j}>=i)=reassign2(outRange2{j}>=i);
    end
end


%%

for i=1:length(ccell)
    idAll=[];
    for j=1:length(outRange2{i})
        cellIdx=outRange2{i}(j);
        idAll=cat(2,idAll,ccell{cellIdx}{outRange2{cellIdx}==i});
        
        
    end
    ccell3{i}=(idAll);
    
end




%% test MCL DIDNT LIKE IT
Tsteady=transitionMatrix;
r=5;
for i=1:100;
    Tsteady=Tsteady*transitionMatrix;
    Tsteady=Tsteady.^r;
    Tsteady=bsxfun(@rdivide, Tsteady,sum(Tsteady));
    Tsteady(isnan(Tsteady))=0;
    plot(sum(Tsteady'));
    drawnow
    pause(.1)
end
%%
tcorr=corr(Tsteady);
cg = clustergram(tcorr);
%%
%tcorr2=corr(transitionMatrix);
tcorr2=corr(transitionMatrix);

tcorr2good=(sum(tcorr2>0))>1;
cg2=clustergram(tcorr2);
cgIdx=str2double(get(cg2,'RowLabels'));
cgIdx2=intersect(cgIdx,find(tcorr2good));
%%
tcorr2=corr(transitionMatrix);
corrmat=tcorr2-eye(size(tcorr2));
corrmat=squareform(corrmat);
Z=linkage(1-corrmat,'complete');
c=cluster(Z,'cutoff',.9,'criterion','distance');
hist(c,1:max(c))

caccum=accumarray(c,ones(size(c)));
caccum=find(caccum<4);
c(ismember(c,caccum))=nan;

%%
for i=1:length(TrackMatrix);
    TrackDatai=TrackData{i};
    TrackDatai=[TrackDatai c(indexAdd(i)+1:indexAdd(i+1))];
    TrackDatai(isnan(TrackDatai(:,end)),:)=[];
    
    TrackData2{i}=TrackDatai(:,[1 2 3 end]);
end


%%
T3=transitionMatrix^2;%+transitionMatrix;
for i=1:1
    plot(transitionMatrix(:,iCheck));
    hold on
    plot(T3(:,iCheck));
    hold off
    ylim([0 .05])
    pause(1)
    for j=1:(length(indexAdd)-1)
        renormgroup= sum(T3((1+indexAdd(j):indexAdd(j+1)),:));
        T3temp= bsxfun(@rdivide,T3((1+indexAdd(j):indexAdd(j+1)),:),renormgroup);
        T3temp(isnan(T3temp))=0;
        T3((1+indexAdd(j):indexAdd(j+1)),:)=T3temp;
        
    end
    % T3(T3<.5)=0;
    T3=bsxfun(@rdivide, T3,sum(T3));
    T3(isnan(T3))=0;
    % T3=T3*transitionMatrix;
    
end



%%

[V,D]=eig(transitionMatrix);
V=real(V);D=real(D);
tSteadState=V*(D>.99)*V';

%%

[Transformed_M, multilevel_ctrl_pts, multilevel_param] = ...
    gmmreg_L2_multilevel(controlmoved, fiducialsAll(10).model(iIdx,:), 3, [4, 0.2, 0.01], ...
    [0.0000008, 0.0000008, 0.0000008],[(1:10)' (1:10)'], 1,1);

%%

[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(T1, T2, 3, [4, 0.2, 0.01], [0.0000008, 0.0000008, 0.0000008], [0 0 0], 1,0);

[Transformed_M3, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(T1, T3, 3, [4, 0.2, 0.01], [0.0000008, 0.0000008, 0.0000008], [0 0 0], 1,0);
[model, multilevel_ctrl_pts, multilevel_param] = ...
    gmmreg_L2_multilevel(model, scene, level, scales,...
    lambdas, fiducial_indices, spring_constant,showFlag);






