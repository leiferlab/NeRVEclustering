function clusterWormTrackCompiler(filePath,fileOutput)
pwd
display(filePath);
display(fileOutput);
imFolder=fileparts(filePath);
if nargin==1
    fileOutput=filePath;
end
%%
dataFolder=(imFolder);
%%

options.thresh1=0.05;
options.minObjSize=50;
options.filterSize=[10,10,4];
    options.method='invdist';
    options.radius=20;
    options.power=1;

param.dim=3;
param.excessive=4;
 param.quiet=1;
 param.difficult=2.e4;
windowSearch=5;

%%
display(imFolder)
[~,pointStats]=compileTrackMatrix(imFolder);

 presentIdx=cellfun(@(x) ~isempty(x),{pointStats.TrackMatrixi},'uniform',0);
presentIdx=find(cell2mat(presentIdx));
N=max(presentIdx);

presentIdx=1:N;
pointStats=pointStats(1:N);




%% do matching of tracks by clustering similarities in matching matrices
TrackMatrix={pointStats.TrackMatrixi};
TrackMatrix=TrackMatrix(1:N);
nObjs=cell2mat(cellfun(@(x) size(x,1),TrackMatrix,'uniform',0)');
indexAdd=[0; cumsum(nObjs(1:N))]';


 
outRange=1:N;%max(1,i-windowSearch):min(length(TrackData),i+windowSearch);
    subTrackMatrix=TrackMatrix(outRange);
    transitionMatrixSize=cellfun(@(x) size(x,1),subTrackMatrix,'uniform',0);
    transitionMatrixSize=sum(cell2mat(transitionMatrixSize));
    
    
%transitionMatrixi=speye(transitionMatrixSize);
%%
indexAdd2=0:100:100*N;
for i=1:N%:length(TrackMatrix)-windowSearch
    %%
    %transitionMatrix=false(size(cell2mat(TrackMatrix(outRange2{1})),1));
%take individual match matrix
        TrackMatrixTemp=pointStats(i).TrackMatrixi;
       % DMatrixTemp=pointStats(i).DMatrixi_x.^2+pointStats(i).DMatrixi_y.^2 ...
        %+pointStats(i).DMatrixi_z.^2;

        if size(TrackMatrixTemp,1)>10 && any(TrackMatrixTemp(:)) 
TrackMatrixTemp=TrackMatrixTemp(:,:);
          validPoints=TrackMatrixTemp(:)>0;
       % [outRangeoverlap,ia, ib]=intersect(outRange2{j},outRange2{i});
      startIdx=indexAdd(outRange(i));
        TrackMatrixTemp=bsxfun(@plus,TrackMatrixTemp , indexAdd2(1:size(TrackMatrixTemp,2)));
      startIdx=indexAdd(outRange(i));
        startPos=startIdx+(1:size(TrackMatrixTemp,1));
        startPos=repmat(startPos,size(TrackMatrixTemp,2),1)';
        transitionIdxX{i}=startPos(validPoints);
        transitionIdxY{i}=TrackMatrixTemp(validPoints);
       % transitionW{i}=DMatrixTemp(validPoints);
%         transitionIdx=sub2ind([transitionMatrixSize,transitionMatrixSize],startPos(validPoints),TrackMatrixTemp(validPoints));
%         transitionIdx=transitionIdx(~isnan(transitionIdx));
%         transitionMatrixIdx{i}=transitionIdx;
end 
end
%%
transitionIdxX=cell2mat(transitionIdxX');
transitionIdxY=cell2mat(transitionIdxY');

nTransitions=length(transitionIdxY);
transitionMatrixi=sparse(transitionIdxX,transitionIdxY,ones(1,nTransitions),...
    transitionMatrixSize,max(indexAdd2),nTransitions);

transitionMatrixi=transitionMatrixi(:,any(transitionMatrixi));
    
%% build training set on frist 
%    transitionMatrixi=double(transitionMatrixi);
%  tic; tcorr3=corrcoef(transitionMatrixi');toc
%      transitionMatrixi=or(transitionMatrixi,speye(size(transitionMatrixi)));
%     transitionMatrixi=double(transitionMatrixi);
%     
nSelectRange=[];
NTrainingRange=min(400,N-1);
nTraining=min(NTrainingRange,N-1);
nSelect=round(2:NTrainingRange/nTraining:NTrainingRange);
for i=1:length(nSelect)-1
nSelectRange{i}=indexAdd(nSelect(i)):indexAdd(nSelect(i)+1);
end
nSelectAdd=cellfun(@(x) length(x), nSelectRange);
nSelectRange=cell2mat(nSelectRange);
%% correlation and cluster, can take up to 10 minutes
subTranstionMatrix=transitionMatrixi(nSelectRange,:);
      tic; tcorr2=sparseTransitionCorr(subTranstionMatrix',[],1);toc
 %     tic; tcorr2=sparseTransitionCorr(subTranstionMatrix');toc
tsize=size(tcorr2);
 [x,y,v]=find(tcorr2);
 
uTri=x>y;
x=x(uTri);
y=y(uTri);
v=v(uTri);

tcorr2=sparse([x; y], [y;x],[v;v],tsize(1),tsize(2));


% transitionMatrixi=normr(transitionMatrixi);
% subTranstionMatrix=transitionMatrixi(nSelectRange,:);

normTransitionMatrixi=bsxfun(@rdivide,transitionMatrixi,sqrt(sum(transitionMatrixi.^2,2)));
normTransitionMatrixi(isnan(normTransitionMatrixi))=0;
     %% calculate eigenvalues, approximately 100 eigenvalues / minute
     if 0
    tic
    opts.issym=1;
    opts.isreal=1;
    opts.disp=0;
    opts.maxit=500;
    
    opts.v0=(tcorr2(:,1));
    [V,D,flag]=eigs(tcorr2,800,'la',opts);flag
    toc     
    D=sum(D,1);
    Dplot=cumsum(D)';
    dD=gradient(Dplot,5);
    ddD=gradient(dD,5);
    KD=abs(ddD)./(1+dD.^2).^(3/2);
     end
%% 
%tcorr2=full(tcorr2);
   % tcorr2=tcorr2.*~speye(size(tcorr2));
    
   %tcorr2(speye(size(tcorr2))>0)=0;
   tic
   tcorr2_lin=tcorr2(tril(true(length(tcorr2)),-1));
  % tcorr2_lin=tcorr2_lin(:);
   toc;tic
   %tcorr2=squareform(tcorr2)
   Z=linkage(1-tcorr2_lin','complete');
   toc;
    %% cluster
    c=cluster(Z,'cutoff',.65,'criterion','distance'); %normally .9999
    %% raname clusters based on size 
    c=c+1;
    caccum=accumarray(c,ones(size(c))); %how many in each cluster

   caccumN=find(caccum>nTraining*1.2); %bad clusters
    
    c(ismember(c,caccumN))=1;
%the rank of each cluster, giving things a new index based on rank rather than cluster group
[~,iaAccum]=sort(caccum,'descend'); 
    [~,iaAccum]=sort(iaAccum);
    [cUnique,~,ib]=unique(c); % ib is posotion of unique terms)
    cUnique=iaAccum(cUnique);
%     [~,ic]=sort(ia);
%     [~, id]=sort(ic);
    c2=cUnique(ib);
  c2(c==1)=nan;
   c2=c2-1;
    
    

  %%  correction factor, optional, add groups that are cloes but were not clustered
    subTcorr=[];subTcorr2=[];
    uniqueIDs=unique(c2(~isnan(c2)));
    uniqueIDs(uniqueIDs==0)=[];
% average correlations matrices from each subgoup
    for i=1:length(uniqueIDs)
       subTcorr(:,uniqueIDs(i))=mean(tcorr2(:,c2==uniqueIDs(i)),2);
    end
    
        for i=1:length(uniqueIDs)
       subTcorr2(uniqueIDs(i),:)=mean(subTcorr(c2==uniqueIDs(i),:),1);
        end
        
      
        subTcorr2=subTcorr2.*~eye(length(subTcorr2));
        [x,y]=find(triu(subTcorr2)>.35);
        
        for i=1:length(x)
        newSum=sum((c2==x(i)|c2==y(i)));
        if newSum<nTraining*1.2
        c2(c2==x(i))=y(i);
        end
        end

        %% order clusters, remove bad ones. 
                  c2=c2+1;

                caccum=hist(c2,1:max(c2)); %how many in each cluster

                % remove clusteres smaller than 40% of training, or larger
                % than 120%
    caccumN=find(caccum<nTraining*.2|caccum>nTraining*1.2); %bad clusters
    
    c2(ismember(c2,caccumN))=1;
    c2(isnan(c2))=1;
    %% again, raname based on size
    caccum=accumarray(c2,ones(size(c))); %how many in each cluster

                [~,iaAccum]=sort(caccum,'descend'); %the rank of each cluster
    [~,iaAccum]=sort(iaAccum);
    
    [cUnique,ia,ib]=unique(c2); % ib is posotion of unique terms)
    cUnique=iaAccum(cUnique);
    
           c2=cUnique(ib);
  c2(c2==1)=nan;
   c2=c2-1;
     
        
        %%

ccell=mat2cell(c2,nSelectAdd);
    neuronsinFrame=cell2mat(cellfun(@(x) sum(~isnan(x)),ccell,'uniform',0));
    
   %% try to plot the clusters
    close all
    assignedNodes1=find(~isnan(c2));
    c3=c2(~isnan(c2));
    [~,ia]=sort(c3);
    assignedNodes=assignedNodes1(ia);

    %%
    uniqueIDs=unique(c2(~isnan(c2)));
    uniqueIDs(uniqueIDs==0)=[];
    
      masterVec=[];
      masterVecVar=[];
      % find "basis" by averaging over a cluster. 
    for i=1:length(uniqueIDs)
       masterVec(:,uniqueIDs(i))=mean(subTranstionMatrix(c2==uniqueIDs(i),:),1);
       totalVec(:,uniqueIDs(i))=sum(subTranstionMatrix(c2==uniqueIDs(i),:),1);
       masterVecVar(:,uniqueIDs(i))=std(subTranstionMatrix(c2==uniqueIDs(i),:),[],1);
       
    end
    masterVec=bsxfun(@minus,masterVec,mean(masterVec,1));
    
    c3=c2(~isnan(c2));
    caccum=accumarray(c3,ones(size(c3))); %how many in each cluster
    masterVecVar(masterVecVar<.1)=.1;
    masterWeights=1./masterVecVar;
    masterVec=sparse(masterVec);
   % masterVecVar=sparese(masterVecVar);

    masterVec=normc(masterVec.*masterWeights);
    input=subTranstionMatrix;
    input=bsxfun(@rdivide,input,sqrt(sum(input.^2,2)));
    input(isnan(input))=0;
    %masterVec(masterVec<0)=5*masterVec(masterVec<0);
aTest=masterVec'*input';
%aTest(aTest<0)=0;
hitIdx=sub2ind(size(aTest),c3(ia),assignedNodes);
hitmap=false(size(aTest));
hitMap(hitIdx)=true;
hitValues=aTest(hitMap);
notHitValues=aTest(~hitMap);
histX=0:.02:1;
nHit=hist(hitValues,histX);
nMiss=hist(notHitValues,histX);
gamma=nHit./(nHit+nMiss);
hitCutoff=histX(find(gamma>.2,1,'first'));
%hitCutoff=quantile(hitValues,.01);

%% project data onto clustered components and detect values larger than threshold
matchProjections=masterVec'*normTransitionMatrixi';
% %%
% acorrTest=aTest*aTest';
% [V,D]=eig(full(acorrTest));
% aTestV=V'*aTest;
% %%
matchProjections=matchProjections.*(matchProjections>hitCutoff);
%normMatchProjections=bsxfun(@rdivide, matchProjections,sum(matchProjections));
%normMatchProjections(isnan(normMatchProjections))=0;

aMax = max(matchProjections);
aMax(aMax==0)=10;
matchBool=bsxfun(@eq, matchProjections,aMax);
matchBool=matchProjections.*matchBool;
[c4,c4x,c4w]=find(matchBool);
c4=sparse(c4x,ones(size(c4x)),c4,max(indexAdd),1);
c4=full(c4);
c4(c4==0)=nan;
ccell=mat2cell(c4,diff(indexAdd));
ccellProj=sparse(c4x,ones(size(c4x)),matchProjections(matchBool>0),max(indexAdd),1);

ccellProj=mat2cell(full(ccellProj),diff(indexAdd));


%% clear doubles
ccellDoubles=cellfun(@(x) check4doubles(x),ccell,'uniform',0);
frames2check=find(cellfun(@(x) ~isempty(x),ccellDoubles));

for iCheck=(frames2check)'
    ccellTemp=ccell{iCheck};
    doubleTerm=ccellDoubles{iCheck};
    weights=ccellProj{iCheck};
    for iiCheck=doubleTerm'
        doubleSearch=find(ccellTemp==iiCheck);
        w=weights(doubleSearch);
        idx2Remove=doubleSearch(w~=max(w));
        ccellTemp(idx2Remove)=nan;
        weights(idx2Remove)=0;
    end
    ccell{iCheck}=ccellTemp;
    ccellProj{iCheck}=weights;
end

%     masterVec=sparse(masterVec);
%     masterVec=bsxfun(@rdivide,masterVec,sum(masterVec));
%     input=bsxfun(@rdivide,subTranstionMatrix,sum(subTranstionMatrix,2));
    
    
%% add new identities into track matrix, kill off unidentified cell

pointStats2=pointStats;

for i=1:length(ccell)
    
   pointStats2(presentIdx(i)).trackIdx=ccell{i};
   pointStats2(presentIdx(i)).trackWeights=ccellProj{i};
    
    
end
pointStats=pointStats2;
%% YOU SHOULD SAVE HERE %%

save([fileOutput],'pointStats2');
