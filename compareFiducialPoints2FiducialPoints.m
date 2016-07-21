%%%% compare two sets of fiducial point cells by finding correspondences in
%%%% all frames
display('Select both fiducial point files')
fFiles=uipickfiles();
%%
fiducials1=load(fFiles{1});
fiducials2=load(fFiles{2});
fiducials1=fiducials1.fiducialPoints;
fiducials2=fiducials2.fiducialPoints;


%% tracking params
param.dim=3;
param.excessive=4;
 param.quiet=1;
 param.difficult=5.e3;
 param.good=2;
 
 minDist=10;

%%
searchList=find(cellfun(@(x) ~isempty(x),fiducials1));
searchList=searchList(searchList<length(fiducials2));
matchMatrix=nan(150,max(searchList));
compiledTracks=cell(1,length(searchList));
distanceAll=compiledTracks;
fiducials1List=fiducials1(searchList);
fiducials2List=fiducials2(searchList);
parfor iFrame=1:length(searchList)
    display(['starting round' num2str(iFrame)])
    %%
    points1=fiducials1List{iFrame};
    points2=fiducials2List{iFrame};
    if ~isempty(points1)
        empty1=cellfun(@(x) isempty(x),points1);
        empty1list=find(all(~empty1,2));
        points1(empty1)={nan};
        empty2=cellfun(@(x) isempty(x),points2);
        empty2list=find(all(~empty2,2));
        points2(empty2)={nan};  
        points1=cell2mat(points1(:,1:4));
        points2=cell2mat(points2(:,1:4));
        points1=points1(empty1list,:);
        points2=points2(empty2list,:);
        trackIn1=[points1 empty1list ones(size(empty1list))];
        trackIn2=[points2 empty2list 2*ones(size(empty2list))];
        
        trackInput=[trackIn1;trackIn2];
        trackInput=trackInput(:,[1 2 4 5 end]);
        trackInput(:,3)=trackInput(:,3)-min(trackInput(:,3))+1;
        trackInput(:,3)=trackInput(:,3)*10;
        trackOutput=nan;
        trackCounter=0;
        try
        while isnan(trackOutput)
        trackOutput=trackJN(trackInput,20-trackCounter,param);
        trackCounter=trackCounter+1;
        end
         trackPairs=trackOutput(:,4);
        paired1=trackPairs(1:2:end);
         paired2=trackPairs(2:2:end);
         
         matchDistances=trackOutput(1:2:end,1:3)-trackOutput(2:2:end,1:3);
         matchDistanceZ=matchDistances(:,end);

         matchDistances=sqrt(sum(matchDistances.^2,2));
        compiledTracks{iFrame}=[paired1 paired2 matchDistances matchDistanceZ];

        [D,I]=pdist2(points2,points1,'euc','Smallest',1);
        distancePairs=[empty1list,empty2list(I'),D'];
        distanceAll{iFrame}=distancePairs;
        catch me
            me
        end
    end
end


%%
matchDistancesMat=zeros(size(matchMatrix));
matchMatrix2=zeros(size(matchMatrix));
matchDistancesMat2=matchDistancesMat;
for iFrame=1:length(searchList)
    trackPairs=compiledTracks{iFrame};
    distancePairs=distanceAll{iFrame};
    if ~isempty(trackPairs)
            paired1=trackPairs(:,1);
         paired2=trackPairs(:,2);
         matchDistances=trackPairs(:,3);
         matchDistancesZ=trackPairs(:,4);
         dpaired1=distancePairs(:,1);
         dpaired2=distancePairs(:,2);
         ddistance=distancePairs(:,3);
         matchMatrix(paired1,searchList(iFrame))=paired2;
         matchDistancesMat(paired1,searchList(iFrame))=matchDistances;
         matchMatrix2(dpaired1,searchList(iFrame))=dpaired2;
         matchDistancesMat2(dpaired1,searchList(iFrame))=ddistance;
    end
end


        
        
%%
matchMatrixRev=nan(size(matchMatrix));

for i=1:size(matchMatrix,2);
    paired1=find(~isnan(matchMatrix(:,i)));
        if ~isempty(paired1)

    paired2=matchMatrix(paired1,i);
   matchMatrixRev(paired2,i)=paired1;
    end
    
end

%%
dThresh=20;
nMatches=(sum(~isnan(matchMatrix),2)/length(matchMatrix));
matchedCells=nMatches>.5;
matchedCells=find(matchedCells);
matchedTimes=(sum(~isnan(matchMatrix),1))>20;
currentMatrix=matchMatrix(matchedCells,matchedTimes);
currentMatrix2=matchMatrix2(matchedCells,matchedTimes);

subDistMat=matchDistancesMat(matchedCells,matchedTimes);
subDistMat2=matchDistancesMat2(matchedCells,matchedTimes);
neuronLikelyID=mode(currentMatrix,2);
corrMatch=bsxfun(@eq,currentMatrix,neuronLikelyID);
corrMatch2=bsxfun(@eq,currentMatrix2,neuronLikelyID);
subDistMat(corrMatch2 & ~corrMatch)=subDistMat2(corrMatch2 & ~corrMatch);
corrMatch=(corrMatch2|corrMatch) & subDistMat<dThresh;
[~,ia]=sort(sum(corrMatch,2),'descend');

wrongMatch=~corrMatch & ~isnan(currentMatrix) & subDistMat<dThresh;
viewMatch=corrMatch*2+wrongMatch;
viewMatchSorted=sort(viewMatch,2);
viewMatchSorted=viewMatchSorted(ia,:);
imagesc(viewMatchSorted);

%%
matchIDs=mode(matchMatrix,2);
%% lets go back and calculate distances to actual matches
for iFrame=1:length(searchList)
    
      display(['starting round' num2str(iFrame)])
    %%
    points1=fiducials1List{iFrame};
    points2=fiducials2List{iFrame};
    if ~isempty(points1)
        empty1=cellfun(@(x) isempty(x),points1);
        empty1list=find(all(~empty1,2));
        points1(empty1)={nan};
        empty2=cellfun(@(x) isempty(x),points2);
        empty2list=find(all(~empty2,2));
        points2(empty2)={nan};  
        points1=cell2mat(points1(:,1:4));
        points2=cell2mat(points2(:,1:4));
         trackPairs=compiledTracks{iFrame};
         cutoffLength=length(points2);
           if ~isempty(trackPairs) 
paired1=trackPairs(:,1);
paired2=matchIDs(paired1);
paired1(paired2>cutoffLength)=[];
paired2(paired2>cutoffLength)=[];
        trueDistance=points2(paired2,1:3)-points1(paired1,1:3);
        trueDistance=sqrt(sum(trueDistance.^2,2));
        trueDistancesMat(paired1,searchList(iFrame))=trueDistance;
           end
    end
end



%%
subTrueMat=trueDistancesMat(matchedCells,matchedTimes);
wrongD=subDistMat(wrongMatch);
    wrongTrue=subTrueMat(wrongMatch);
    
    %%
    wrongMatch=~corrMatch & ~isnan(currentMatrix) & subDistMat<dThresh;
    wrongMatch=wrongMatch & subTrueMat>subDistMat;
    corrMatch2=(subTrueMat<subDistMat) |corrMatch;
viewMatch=corrMatch2*2+wrongMatch;
%viewMatch(subTrueMat<subDistMat)=3;
viewMatchSorted=sort(viewMatch,2,'descend');
[~,ia]=sort(sum(viewMatch>=2,2),'descend');

viewMatchSorted=viewMatchSorted(ia,:);
imagesc(viewMatchSorted);

%%
cmap=parula;%pmkmp(64,'Swtth');
sampleVolume=fiducials1{100};
samplePoints=cell2mat(sampleVolume);
accuracy=nan(1,length(samplePoints));
accuracy(matchedCells)=1-mean(corrMatch2,2);
accuracy(accuracy>.15)=.15;
% accuracy=normalizeRange(accuracy);
colorAccuracy=squeeze(ind2rgb(round(accuracy*62+2),cmap));
colorAccuracy(isnan(accuracy),:)=.9;

scatter3sph(samplePoints(:,1),samplePoints(:,2),samplePoints(:,3)*30,...
    'size',5,'color',squeeze(colorAccuracy));axis equal off
colorbar;colormap(cmap);caxis([0 .15]);