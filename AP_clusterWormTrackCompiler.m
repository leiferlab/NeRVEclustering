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
PS_ref=load([dataFolder filesep 'pointStatsRef']);
PS_ref=PS_ref.PS_ref;
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
n_ref_neurons=cellfun(@(x) size(x,1),{PS_ref.straightPoints});

%% loop through all track matrices to build binary feature vectors
%each track matrix has an n neurons by t time points vector with each
%neuron having a number (match index) for each time. To build the binary
%feature vector 

%build binary sparce matrix, referred to as a transition matrix
transitionMatrixCell=...
    cellfun(@(x) oneHotNeuron(x,n_ref_neurons),TrackMatrix,'Uniform',0);

%% build training set on first 800 time points
%if this fails, continuously decrease the number of training time points
%until it works. 
% start later in the recording -- in case there are 'bad volumes'
nStart= 500
for iTry=0:5
    %% try loop to try to work around out of memory issues

        NTrainingRange=900-100*iTry+nStart;
	%bla
	%NTrainingRange=450-50*iTry;	
        
        NTrainingRange=min(NTrainingRange,N-1);
        nTraining=min(NTrainingRange-nStart,N-1-nStart);
        
        display(['Attempting nTraining of ' num2str(nTraining)]);
        nSelect=round(nStart:NTrainingRange/nTraining:NTrainingRange);
    try
        
        %% correlation and cluster, can take up to 10 minutes  
        %make subset of transition matrix
        subTranstionMatrix=cell2mat(transitionMatrixCell(nSelect)');
        %calculate distance matrix as correlation distance using custom 
        %code for sparce correlations.
        tic;
        tcorr2=sparseTransitionCorr(subTranstionMatrix',[],1);
        display(['Distance calculated in ' num2str(toc) 's']);
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

%%bla
n_clusters

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

%bla
size(caccumN)

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

%%bla
cluster_items

cluster_items(cluster_items==0)=[];
n_clusters=max(cluster_items);

%%bla
n_clusters

display(['Remaining number of clusters ' num2str(n_clusters)]);
subTcorr=zeros(length(tcorr2),n_clusters);




%TODO: trying this
subTcorr2=zeros(n_clusters.').';
%subTcorr2=zeros(n_clusters);



% average correlations matrices from each subgoup
for i=1:n_clusters
    i_cluster=cluster_items(i);
    subTcorr(:,i_cluster)=mean(tcorr2(:,cluster_assign==i_cluster),2);
end

for i=1:n_clusters
    i_cluster=cluster_items(i);
    subTcorr2(i_cluster,:)=mean(subTcorr(cluster_assign==i_cluster,:),1);
end


%% find "basis" by averaging over a cluster. this is the cm of the point clound
% and will be used to classify every neuron from every volume, this is step
% 2 of the clustering in the paper. 

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

%make list of indices that should be theoretically assigned to the group
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
output=cellfun(@(x) distanceClassify(masterVec,x,hitCutoff),...
    transitionMatrixCell','Uniform',0);
output=cat(1,output{:});
track_cell=output(:,1);
weight_cell=output(:,2);
%% add new identities into track matrix, kill off unidentified cell

pointStats2=pointStats;
for i=1:length(pointStats2)
    pointStats2(presentIdx(i)).trackIdx=track_cell{i};
    pointStats2(presentIdx(i)).trackWeights=weight_cell{i};
end
%% SAVE HERE %%
fileOutput_stats=strrep(fileOutput,'.mat','_info.mat');
save(fileOutput,'pointStats2');
save(fileOutput_stats,'masterVec','hitCutoff'...
    ,'cluster_assign','caccum','subTcorr2')


%% write time end stamp
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
if contains(hostname,'della')
    Fid=fopen([dataFolder filesep 'status.txt'],'a');
    status=[datestr(datetime('now')) ': Finished NERVE tracking \n'];
    fprintf(Fid,status);
    fclose(Fid);
end
