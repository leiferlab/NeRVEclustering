function clusterWormTrackCompiler(filePath,fileOutput)

% Run after all images have straightened and segmented, and all the track
% matrices have been created. This program is built to run on della. It
% takes all of the positional fingerprints (track matrices) and clusters
% them. The input is the file path to the main folder, which has a folder
% with all of the track matrices (from clusterWormTracker). It also
% requires a pointsstats.mat file in the folder. The program has no output,
% but ssaves a pointStats2.mat folder, which adds a trackIdx field to the
% original pointStats.

% get folder
if nargin==1
    fileOutput=filePath;
end
display(filePath);
display(fileOutput);


%%
dataFolder=fileparts(filePath);

if isempty(dataFolder)
    dataFolder=pwd;
end

display(['Parent Folder is ' dataFolder]);
display([ 'PS file is ' filePath]);


%load all pointStats with, which has all trackMatrix data
[~,pointStats]=compileTrackMatrix(dataFolder);

% fill in indexing holes with error fits,
for iPS=1:length(pointStats)
    if isempty(pointStats(iPS).stackIdx);
        pointStats(iPS).stackIdx=iPS;
    end
end

%find which pointStats had tracking completed
track_present=cellfun(@(x) ~isempty(x),{pointStats.TrackMatrixi},'uniform',0);
N=find(cell2mat(track_present),1,'last');
display([ 'Number of frames is ' num2str(N)]);

%only go up to last trackmatrix present
pointStats=pointStats(1:N);
presentIdx=1:length(pointStats);




%% do matching of tracks by clustering similarities in matching matrices
%build cell array of all track matrices
TrackMatrix={pointStats(presentIdx).TrackMatrixi};
%get number of neurons found in each volume
nObjs=cell2mat(cellfun(@(x) size(x,1),TrackMatrix,'uniform',0)');

%neurons will all be linearly indexed, so a number is added to each index
%depending on the timepoint the neuron is found in.
indexAdd=[0; cumsum(nObjs(1:N))]';
%for the references this can be pre
indexAdd2=0:100:100*N;

transitionMatrixSize=sum(nObjs);
%% loop through all track matrices to build binary feature vectors
%each track matrix has an n neurons by t time points vector with each
%neuron having a number (match index) for each time. To build the binary
%feature vector 
    transitionIdxX=cell(1,N);
    transitionIdxY=cell(1,N);
    
for i=1:N
    %take individual match matrix and linearize all matches
    TrackMatrixTemp=pointStats(i).TrackMatrixi;

    if size(TrackMatrixTemp,1)>10 && any(TrackMatrixTemp(:))
        validPoints=TrackMatrixTemp(:)>0;
        %add index offset for for reference time points
        indexAdd_i= indexAdd2(1:size(TrackMatrixTemp,2));
        TrackMatrixTemp=bsxfun(...
            @plus,TrackMatrixTemp ,indexAdd_i);
        %add indexing offset for sample
        startIdx=indexAdd(presentIdx(i));
        startPos=startIdx+(1:size(TrackMatrixTemp,1));
        startPos=repmat(startPos,size(TrackMatrixTemp,2),1)';
        %get indeces for binary feature vectors
        transitionIdxX{i}=startPos(validPoints);
        transitionIdxY{i}=TrackMatrixTemp(validPoints);
    end
end
%% 
%convert results into matrices
transitionIdxX=cell2mat(transitionIdxX');
transitionIdxY=cell2mat(transitionIdxY');

%build binary sparce matrix, referred to as a transition matrix
nTransitions=length(transitionIdxY);
transitionMatrixi=sparse(transitionIdxX,transitionIdxY,ones(1,nTransitions),...
    transitionMatrixSize,max(indexAdd2),nTransitions);

%collapse away empty columns
transitionMatrixi=transitionMatrixi(:,any(transitionMatrixi));

%% build training set on first 800 time points
%if this fails, continuously decrease the number of training time points
%until it works. 

for iTry=0:5
    %% try loop to try to work around out of memory issues

        nSelectRangeCell=[];
        NTrainingRange=800-100*iTry;
        NTrainingRange=min(NTrainingRange,N-1);
        nTraining=min(NTrainingRange,N-1);
        
        display(['Attempting nTraining of ' num2str(nTraining)]);
        nSelect=round(1:NTrainingRange/nTraining:NTrainingRange);
        nSelectRangeCell=cell(1,length(nSelect)-1);
        for i=1:length(nSelect)-1
            nSelectRangeCell{i}=1+indexAdd(nSelect(i)):indexAdd(nSelect(i)+1);
        end
        nSelectAdd=cellfun(@(x) length(x), nSelectRangeCell);
        nSelectRange=cell2mat(nSelectRangeCell);
        
    try
        
        %% correlation and cluster, can take up to 10 minutes  
        %make subset of transition matrix
        subTranstionMatrix=transitionMatrixi(nSelectRange,:);
        %calculate distance matrix as correlation distance using custom 
        %code for sparce correlations.
        tic;
        tcorr2=sparseTransitionCorr(subTranstionMatrix',[],1);
        display(['Distance calculated in ' num2str(toc) 's']);
        transition_mag=sqrt(sum(transitionMatrixi.^2,2));
        normTransitionMatrixi=bsxfun(...
            @rdivide,transitionMatrixi,transition_mag);
        normTransitionMatrixi(isnan(normTransitionMatrixi))=0;
   
        %get lower triangular section of the matrix for using clustering
        tcorr2_lin=tcorr2(tril(true(length(tcorr2)),-1));
        
        %% cluster
        %find clusters with a cutoff distance .9
        Z=linkage(1-tcorr2_lin','complete');
        %clusters_initial is a linear list for each neuron in the subtransition matrix
        %with a number indicating which cluster that neuron belongs to
        clusters_initial=cluster(Z,'cutoff',.9,'criterion','distance'); 
        
        %if clusters successfully calculated, exit loop. 
        display(['Success with nTraining of ' num2str(nTraining)]);
        break
    catch me
        display(['Failed with nTraining of ' num2str(nTraining)]);
        me
    end
    
end

%% remove nans and zeros
cluster_items=unique(clusters_initial(~isnan(clusters_initial)));
cluster_items(cluster_items==0)=[];
n_clusters=max(cluster_items);

%%  correction factor, add groups that are close but were not clustered

%build a new distance matrix with distance between formed clusters. 
subTcorr=zeros(length(tcorr2),n_clusters);
subTcorr2=zeros(n_clusters);
% average correlations matrices from each subgoup in each of the two
% directions
for i=1:n_clusters
    i_cluster=cluster_items(i);
    subTcorr(:,i_cluster)=mean(tcorr2(:,clusters_initial==i_cluster),2);
end

for i=1:n_clusters
    i_cluster=cluster_items(i);
    subTcorr2(i_cluster,:)=mean(subTcorr(clusters_initial==i_cluster,:),1);
end
%get value dy the diagnol, use it as distribution to determine whether off
%diagnol terms are significant.
inGroup=subTcorr2(eye(n_clusters)>0);
%make a cutoff to determine if diagnol terms should be significant. 
regroupCutoff=max(mean(inGroup)-std(inGroup)*2,.25);

%find clusters matched above cutoff
[x,y,v]=find(subTcorr2.*((triu(subTcorr2,1)>regroupCutoff)));
[~,ia]=sort(v,'descend');
x=x(ia);
y=y(ia);
%% recombine clusters that are determined to have been matched together. 

for i=1:length(x)
    % get the neurontimes that belong into the clusters to be merged
    newSumx=sum(clusters_initial==x(i));
    newSumy=sum(clusters_initial==y(i));
    %the new number in that cluster is the sum of the two old clusters
    newSum=newSumx+newSumy;
    %if the sum of the two clusters is not more than the number of training
    %times x 1.2, just straight merge them
    if newSum<nTraining*1.2
        clusters_initial(clusters_initial==x(i))=y(i);
    else
        %if they are larger combine the two clusters and recluster them so
        %that the number of objects is less then 1.2xnTraining. 
        checkx=find(clusters_initial==x(i));
        checky=find(clusters_initial==y(i));
        check_set=[checkx;checky];
        %make distance matrix for subset
        checkCorr=tcorr2(check_set,check_set);
        checkCorr=checkCorr(tril(true(newSum),-1));
        %make clusters
        Ztest=linkage(1-full(checkCorr)','complete');
        recluster=cluster(Ztest,'maxclust',2);
        if max(accumarray(recluster,ones(1,newSum)))<nTraining*1.2
            clusters_initial(check_set(recluster==1))=x(i);
            clusters_initial(check_set(recluster==2))=y(i);
        end
        
        
    end
end

%% order clusters, remove bad ones.
%add one to clusters, anything equal to one will be removed
clusters_initial=clusters_initial+1;
 %how many in each cluster
caccum=hist(clusters_initial,1:max(clusters_initial));
% remove clusteres smaller than 40% of training, or larger
% than 120%
caccumN=find(caccum<nTraining*.4|caccum>nTraining*1.2); %bad clusters

% set bad clusters and nans to one for later removal. 
clusters_initial(ismember(clusters_initial,caccumN))=1;
clusters_initial(isnan(clusters_initial))=1;
%%  raname based on size of cluster, so cluster 1 will be largest
caccum=accumarray(clusters_initial,ones(size(clusters_initial))); 

[~,iaAccum]=sort(caccum,'descend'); %the rank of each cluster
[~,iaAccum]=sort(iaAccum);
% ib is posotion of unique terms)
%cluster_items are the unique labels for the clusters
[cluster_items,~,ib]=unique(clusters_initial); 
cluster_items=iaAccum(cluster_items);

cluster_assign=cluster_items(ib);
cluster_assign(cluster_assign==1)=nan;
cluster_assign=cluster_assign-1;

%%  make new distance matrix between clusters for saving, 
%not currently used elsewhere
cluster_items=unique(cluster_assign(~isnan(cluster_assign)));
cluster_items(cluster_items==0)=[];
n_clusters=max(cluster_items);

subTcorr=zeros(length(tcorr2),n_clusters);
subTcorr2=zeros(n_clusters);
% average correlations matrices from each subgoup
for i=1:n_clusters
    i_cluster=cluster_items(i);
    subTcorr(:,i_cluster)=mean(tcorr2(:,cluster_assign==i_cluster),2);
end

for i=1:n_clusters
    i_cluster=cluster_items(i);
    subTcorr2(i_cluster,:)=mean(subTcorr(cluster_assign==i_cluster,:),1);
end


%% find "basis" by averaging over a cluster.

masterVec=zeros(n_clusters,size(subTranstionMatrix,2));
totalVec=zeros(n_clusters,size(subTranstionMatrix,2));
masterVecVar=zeros(n_clusters,size(subTranstionMatrix,2));
% 
for i=1:n_clusters
    %loop over clusters get average, sum and std of the trainingset
    %transition matrix for the items in a cluster
    i_cluster=cluster_items(i);
    masterVec(i_cluster,:)=mean(...
        subTranstionMatrix(cluster_assign==i_cluster,:),1);
    totalVec(i_cluster,:)=sum(...
        subTranstionMatrix(cluster_assign==i_cluster,:),1);
    masterVecVar(i_cluster,:)=std(....
        subTranstionMatrix(cluster_assign==i_cluster,:),[],1);
    
end
%subtract off mean to get mean subtracted means of the trainingset
%transition matrix
masterVec=bsxfun(@minus,masterVec,mean(masterVec,1));
%make a minimum value for the std
masterVecVar(masterVecVar<.1)=.1;
%use 1/std as the weighting for each vector
masterWeights=1./masterVecVar;
masterVec=sparse(masterVec);

%apply the weightings to the mastervectors
masterVec=normr(masterVec.*masterWeights);



%% project training set data onto master vectors for each cluster
input=subTranstionMatrix;
input=bsxfun(@rdivide,input,sqrt(sum(input.^2,2)));
input(isnan(input))=0;

%project training data onto masters, score is a n training neurons by m
%clusters, with the match score of each cluster in the rows. 
score=masterVec*input';




%% 
%compare the distribution of values in the projection that are expected to
%be  hits and the ones expected not to be. Find the cutoff which makes us
%80% sure that a hit should be a hit. 

assignedNodes=find(~isnan(cluster_assign));
cluster_assign_present=cluster_assign(~isnan(cluster_assign));
[~,ia]=sort(cluster_assign_present);
assignedNodes=assignedNodes(ia);

%make list of indices that should be theoretically 
hitIdx=sub2ind(size(score),cluster_assign_present(ia),assignedNodes);
%turn that list onto a logical matrix, where each training neuron has 1
%for the cluster it belongs to and 0 everywhere else. 
hitMap=false(size(score));
hitMap(hitIdx)=true;
% get the actual values of the projections at the locations we expect 1
hitValues=score(hitMap);
% and get the other values to compare
notHitValues=score(~hitMap);

%compare distributions to find cutoff
histX=0:.02:1;
%histogram two distributions
nHit=hist(hitValues,histX);
nMiss=hist(notHitValues,histX);
%find probability that a value is a hit
gamma=nHit./(nHit+nMiss);
%cutoff is where p(hit) is greater than 20%. 
hitCutoff=histX(find(gamma>.2,1,'first'));

%% project data onto clustered components and detect values larger than threshold
matchProjectionsRaw=masterVec*normTransitionMatrixi';
matchProjections=matchProjectionsRaw.*(matchProjectionsRaw>hitCutoff);

aMax = max(matchProjections);
%if there is no hit, set the max high, the projection will never equal it
aMax(aMax==0)=10;
%find where the best match occurs
matchBool=bsxfun(@eq, matchProjections,aMax);
matchBool=matchProjections.*matchBool;
[match_x,match_y,~]=find(matchBool);

% make sparse linear matrix with all matches, match_x has the clustered index 
% for each neuron 
match_x=sparse(match_y,ones(size(match_y)),match_x,max(indexAdd),1);
match_x=full(match_x);
match_x(match_x==0)=nan;

%turn all cluster assignments into cell
ccell=mat2cell(match_x,diff(indexAdd));
%turn all cluster weights (value of projection onto that cluster master
%vec) into cell
ccellProj=sparse(...
    match_y,ones(size(match_y)),matchProjections(matchBool>0),max(indexAdd),1);
ccellProj=mat2cell(full(ccellProj),diff(indexAdd));
matchProjectionsCell=mat2cell(...
    full(matchProjectionsRaw),n_clusters,diff(indexAdd));


%% clear doubles
%find if there are time points with doubles
ccellDoubles=cellfun(@(x) check4doubles(x),ccell,'uniform',0);
frames2check=find(cellfun(@(x) ~isempty(x),ccellDoubles));

%loop over all frames with doubles
for iCheck=(frames2check)'
    %get the indices which appear twice
    ccellTemp=ccell{iCheck};
    doubleTerm=ccellDoubles{iCheck};
    %also get the weights for each of them, 
    weights=ccellProj{iCheck};
    %loop over all doubles in that frame
    for iiCheck=doubleTerm'
        % get the neurons that were assigned to the same index
        doubleSearch=find(ccellTemp==iiCheck);
        % get the weights for each of them
        w=weights(doubleSearch);
        % remove the one with the smaller projection. 
        idx2Remove=doubleSearch(w~=max(w));
        ccellTemp(idx2Remove)=nan;
        % also set that weight to zero
        weights(idx2Remove)=0;
    end
    %repopulate the cells
    ccell{iCheck}=ccellTemp;
    ccellProj{iCheck}=weights;
end


%% add new identities into track matrix, kill off unidentified cell

pointStats2=pointStats;
for i=1:length(ccell)
    pointStats2(presentIdx(i)).trackIdx=ccell{i};
    pointStats2(presentIdx(i)).trackWeights=ccellProj{i};
end
%% SAVE HERE %%
fileOutput_stats=strrep(fileOutput,'.mat','_info.mat');
save(fileOutput,'pointStats2');
save(fileOutput_stats,'masterVec','matchProjectionsCell','cInitial'...
    ,'cluster_assign','caccum','subTcorr2')
