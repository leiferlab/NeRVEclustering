

imFolder=uipickfiles%'F:\20141212\BrainScanner20141212_145951\CLstraight_20150525_2';
imFolder=imFolder{1};
%%
zindexer=@(x,s) x./(s)+1;

options.thresh1=0.05;
options.minObjSize=50;
options.filterSize=[10,10,4];
    options.method='invdist';
    options.radius=20;
    options.power=1;
imageRange=500:750;

param.dim=3;
param.excessive=4;
 param.quiet=1;
 param.difficult=2.e4;
windowSearch=5;

%%
[TrackMatrix,pointStats]=compileTrackMatrix2(imFolder);

 presentIdx=cellfun(@(x) ~isempty(x),{pointStats.TrackMatrixi},'uniform',0);
presentIdx=find(cell2mat(presentIdx));
N=max(presentIdx);

presentIdx=1:N;
pointStats=pointStats(1:N);



%% do matching of tracks by clustering similarities in matching matrices
TrackMatrix={pointStats.TrackMatrixi};
TrackMatrix=TrackMatrix(1:N);
indexAdd=cell2mat(cellfun(@(x) size(x,1),TrackMatrix,'uniform',0)');
indexAdd=[0; cumsum(indexAdd(1:N))]';


 
outRange=1:N;%max(1,i-windowSearch):min(length(TrackData),i+windowSearch);
    subTrackMatrix=TrackMatrix(outRange);
    transitionMatrixSize=cellfun(@(x) size(x,1),subTrackMatrix,'uniform',0);
    transitionMatrixSize=sum(cell2mat(transitionMatrixSize));
    
    
%transitionMatrixi=speye(transitionMatrixSize);
transitionMatrixIdx=[];
%%
clear transitionIdxX transitionIdxY
for i=1:N%:length(TrackMatrix)-windowSearch
    %%
    %transitionMatrix=false(size(cell2mat(TrackMatrix(outRange2{1})),1));
%take individual match matrix
        TrackMatrixTemp=TrackMatrix{i};
        if ~isempty(TrackMatrixTemp)

TrackMatrixTemp=TrackMatrixTemp(:,1:N);
                        validPoints=TrackMatrixTemp(:)>0;
       % [outRangeoverlap,ia, ib]=intersect(outRange2{j},outRange2{i});
      startIdx=indexAdd(outRange(i));
        TrackMatrixTemp=bsxfun(@plus,TrackMatrixTemp , indexAdd(1:N));

        startPos=startIdx+(1:size(TrackMatrixTemp,1));
        startPos=repmat(startPos,N,1)';
        transitionIdxX{i}=startPos(validPoints);
        transitionIdxY{i}=TrackMatrixTemp(validPoints);
%         transitionIdx=sub2ind([transitionMatrixSize,transitionMatrixSize],startPos(validPoints),TrackMatrixTemp(validPoints));
%         transitionIdx=transitionIdx(~isnan(transitionIdx));
%         transitionMatrixIdx{i}=transitionIdx;
end 
end
i=0;
transitionIdxX=cell2mat(transitionIdxX');
transitionIdxY=cell2mat(transitionIdxY');
nTransitions=length(transitionIdxY);
transitionMatrixi=sparse(transitionIdxX,transitionIdxY,ones(1,nTransitions),...
    transitionMatrixSize,transitionMatrixSize,nTransitions);

%%


% %%
% for i=1:N
%     transitionMatrixi(transitionMatrixIdx{i})=true;
%     
%     
% end

    %%
  %  transitionMatrixi=transitionMatrixi+transitionMatrixi'+eye(size(transitionMatrixi));
   % transitionMatrixi=bsxfun(@rdivide, transitionMatrixi,sum(transitionMatrixi));
    %transitionMatrixi(isnan(transitionMatrixi))=0;
    %%
%    transitionMatrixi=double(transitionMatrixi);
 %  tic; tcorr3=corrcoef(transitionMatrixi');toc
      tic; tcorr2=sparseTransitionCorr(transitionMatrixi');toc

      %%
%tcorr2=full(tcorr2);
    corrmat=tcorr2.*~speye(size(tcorr2));
    corrmat=squareform(corrmat);
    
    Z=linkage(1-corrmat,'complete');
    
    c=cluster(Z,'cutoff',.9999,'criterion','distance');
    
    caccum=accumarray(c,ones(size(c)));
    caccum=find(caccum<100 | caccum>N*1.2);
    c(ismember(c,caccum))=0;
    [~,ia,ib]=unique(c);
    [~,ic]=sort(ia);
    [~, id]=sort(ic);
    c2=id(ib);
    c2(c==0)=nan;

  %%  correction factor, optional
    subTcorr=[];subTcorr2=[];
    uniqueIDs=unique(c2(~isnan(c2)));
    uniqueIDs(uniqueIDs==0)=[];
    for i=1:length(uniqueIDs)
       subTcorr(:,uniqueIDs(i))=mean(tcorr2(:,c2==uniqueIDs(i)),2);
    end
    
        for i=1:length(uniqueIDs)
       subTcorr2(uniqueIDs(i),:)=mean(subTcorr(c2==uniqueIDs(i),:),1);
        end
    
        subTcorr2=subTcorr2.*~eye(length(subTcorr2));
        [x,y]=find(triu(subTcorr2)>.5);
        
        for i=1:length(x)
            c2(c2==x(i))=y(i);
        end
        
        %%
          
    caccum=hist(c2,1:max(c2));%accumarray(c2,ones(size(c2)));
    caccum=find(caccum<200 | caccum>N*1.2);
    c2(ismember(c,caccum))=nan;
            ccell=mat2cell(c2,diff(indexAdd(min(outRange):max(outRange+1))));
    neuronsinFrame=cell2mat(cellfun(@(x) sum(~isnan(x)),ccell,'uniform',0));
    
    
    
   %%
    
    assignedNodes1=find(~isnan(c2));
    c3=c2(~isnan(c2));
    [~,ia]=sort(c3);
    assignedNodes=assignedNodes1(ia);
    imagesc(1-tcorr2(assignedNodes,assignedNodes))
    caxis([0,1])

%% add new identities into track matrix, kill off unidentified cell

pointStats2=pointStats;

for i=1:length(ccell)
    
   pointStats2(presentIdx(i)).trackIdx=ccell{i};
    
    
    
end
pointStats=pointStats2;

%%
pointStats3=pointStats2;
%for iIteration=1:3
for i =1:N
    TrackMatrix2=zeros(size(pointStats2(i).TrackMatrixi));
    if ~isempty(TrackMatrix2) && ~isempty(pointStats2(i).trackIdx) && neuronsinFrame(i)>30
    for j=1:N
        TrackMatrix2Temp=pointStats2(i).TrackMatrixi(:,j);
        TrackMatrix2Lookup=pointStats2(j).trackIdx;
        trackLoc=find(TrackMatrix2Temp);
        TrackMatrix2Temp2=TrackMatrix2Temp(trackLoc);
        TrackMatrix2Temp2=TrackMatrix2Lookup(TrackMatrix2Temp2);
        TrackMatrix2Temp(trackLoc)=TrackMatrix2Temp2;
       TrackMatrix2(:,j)= TrackMatrix2Temp;

    end
    TrackMatrix2(TrackMatrix2==0)=nan;
    trackMode=mode(TrackMatrix2,2);
    modeNum=mean(bsxfun(@eq,TrackMatrix2,trackMode),2);
    newPoints=modeNum>.2 & isnan(pointStats2(i).trackIdx) & ~ismember(trackMode,pointStats3(i).trackIdx);
    [find(newPoints) pointStats2(i).trackIdx(newPoints) trackMode(newPoints)];
    pointStats3(i).trackIdx(newPoints)=trackMode(newPoints);
    end
end
%end



pointStats2=pointStats3;


%% match control points with automated points


for i=1:length(pointStats2)
   %      pointStats2(presentIdx(i)).transitionMatrix=TrackMatrix{i};

    try
    rawPoints=pointStats2(presentIdx(i)).rawPoints;
    rawPoints=rawPoints(~isnan(pointStats2(presentIdx(i)).trackIdx),:);
    rawPointsId=pointStats2(presentIdx(i)).trackIdx(~isnan(pointStats2(presentIdx(i)).trackIdx));
    controlPoints=pointStats2(presentIdx(i)).controlPoints;
    rawPoints=[rawPoints rawPointsId ones(size(rawPoints,1),1)];
    controlPoints=[controlPoints (1:size(controlPoints,1))' 2*ones(size(controlPoints,1),1)];
    matchOut=trackJN([rawPoints;controlPoints],9,param);
    matchID=matchOut(:,end-2);
    pairedId=[ matchID(and(matchOut(:,end-1)==1, circshift(~diff([1;matchOut(:,end)]),-1))) ...
        matchID(and(matchOut(:,end-1)==2, ~diff([1;matchOut(:,end)])))];
    trackIdx=pointStats2(presentIdx(i)).trackIdx;
    trackIdx(~ismember(trackIdx,pairedId(:,1)))=nan;
    [~,ia]=sort(trackIdx);
    pairedId=sortrows(pairedId);
    matchIdx=nan(size(pointStats2(presentIdx(i)).trackIdx));
    matchIdx(ia(1:size(pairedId,1)))=pairedId(:,2);
    
    pointStats2(presentIdx(i)).matchIdx=matchIdx;
    catch ME
        display(['Error frame:' num2str(presentIdx(i))]);
    end
    
end

%% turn matches into matrix for checking accuracy
nTracks=max(cell2mat(cellfun(@(x) nanmax(x),{pointStats2.trackIdx},'uniform',0)'));
matchTest=zeros(nTracks,N);
Dmat=nan(size(matchTest));
lookupMat=nan(size(matchTest));
for i=1:N
    trackIdx=pointStats2(i).trackIdx;
    matchIdx=pointStats2(i).matchIdx;
    pointIdx=pointStats2(i).pointIdx;
        if nnz(~isnan(trackIdx))>10
                noNans=~isnan(trackIdx+matchIdx);
    trackIdx=trackIdx(noNans);matchIdx=matchIdx(noNans);
    pointIdx=pointIdx(noNans);
    
 rawPoints=pointStats2(i).rawPoints(noNans,:);
 controlPoints=pointStats2(i).controlPoints(matchIdx,:);
 D=sqrt(sum((rawPoints-controlPoints).^2,2));

    matchTest(trackIdx,i)=matchIdx;
    lookupMat(trackIdx,i)=pointIdx;
    Dmat(trackIdx,i)=D;
        end
    
        
end

%%
% Dmat=Dmat(any(matchTest,2),:);
% lookupMat=lookupMat(any(matchTest,2),:);
% matchTest=matchTest(any(matchTest,2),:);
matchTestRaw=matchTest;
matchErrRaw=bsxfun(@ne,matchTest,mode(matchTest,2)).*~isnan(matchTest);

matchTest(:,neuronsinFrame<40)=nan;
Dmat(matchTest==0)=nan;

matchTest(matchTest==0)=nan;
matchErr=bsxfun(@ne,matchTest,mode(matchTest,2)).*~isnan(matchTest);
matchCorr=bsxfun(@eq,matchTest,mode(matchTest,2));

nCorr=hist(Dmat(~matchErr),0:15);
nErr=hist(Dmat(matchErr>0),0:15);


%%
rawLookupmat=bsxfun(@plus, lookupMat,indexAdd(1:end-1));
subAssignedNodes=assignedNodes;%(assignedNodes<5000);
failedPoints=rawLookupmat(matchErr>0);
failedPoints=assignedNodes(ismember(subAssignedNodes,failedPoints));
hitPoints=assignedNodes(~ismember(subAssignedNodes,failedPoints));
% imagesc(1-tcorr2(subAssignedNodes,subAssignedNodes))
% hold on
% scatter(failedPoints,failedPoints,'rx');
% hold off

%% map errors into structure
for i=1:length(TrackMatrix)
    trackIdx=pointStats2(i).trackIdx;
    matchIdx=pointStats2(i).matchIdx;
    pointIdx=pointStats2(i).pointIdx;
    errIdx=find(matchErrRaw(:,i));
    if ~isempty(errIdx);
        pointStats2(i).errIdx=(ismember(trackIdx,errIdx));
    end
    
    
end

%% show error rate and error map
totErr=sum(matchErr,2);
totCorr=sum(matchTest>0,2)-sum(matchErr,2);
totMissing=sum(isnan(matchTest),2);

Y=[totCorr,totErr,totMissing];
[Y,iy]=(sortrows(Y,1));
Y=flipud(Y);
iy=flipud(iy);
bar(Y,'stacked');

matchMap=double(matchCorr);
matchMap(matchErr>0)=2;
matchMap(isnan(matchTest))=3;
xRange=find(any(~isnan(matchTest)));
figure;
imagesc(matchMap(iy,:))
xlim([min(xRange) max(xRange)]);
figure
bar(Y,'stacked')
ylim([0 N]);
xlim([0 70]);
%%
ccell2=ccell;

%loop over all clusters
for i=1:length(ccell2)-1
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
