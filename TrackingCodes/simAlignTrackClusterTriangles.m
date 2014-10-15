% simAlign.m
% 
% Calculate a metric of alignment between two noisy measurements of the 3D 
% position of neuron cell bodies.
%
% The position of neurons in a worm are measured twice, during acquisition A
% and acquisition B. The measurement is noisy, and the worm may be imaged
% in different orientations.
%
% To aling the two acquisitions, it is important to have a metric of the
% similarity of the two sets of neuron poisitions. This will allow for
% functions that try many different orientations while optimizing for
% similarity.
%
% This is a small script that simulates two noisy 3D neural position
% datasets that are already aligned. It then  pairs each neuron from A 
% with a neuron from B. Finally it measures how aligned the sets of neurons 
% are by measuring the sum of the root mean squared displacement between 
% each pair.
%
% Andrew Leifer
% leifer@princeton.edu
% 8 February 2014


%We start by assuming that the user acquired two 3D fluorescent images 
%of neurons of a worm. We will call them acquisition A & B.


%Number of Neurons Identified 
%The number of neurons aren't necessarily the same because it is hard for a
%  user to tell what exactly is a neuron.

N=120; %acquisition A
M=120; %acquisition B
assert(N>=M); %by convention, the code assumes A has more neurons then B

%set dimensions of  simulated 3D volume in microns
xdiam=75;
ydiam=150;
zdiam=75;

%set standard deviation of jitter in neural position between acquisition A & B
pstd=2;
permSize=40;
numIterations=1;

%Set the markers to use when plotting A & B
Bmarker='rs';
Amarker='bo';

%% Simulate position data 

Apos=[rand(N,1).*xdiam, rand(N,1).*ydiam, rand(N,1).*zdiam];

mapping=randperm(N);  %Elements are neuron ID's from A  and their index is the neuorn ID from B
[~,unmapping]=sort(mapping);
Bpos=Apos(mapping,:); %Shuffle the order of neurons from A and assing to B
Bpos=Bpos+randn(size(Bpos)).*pstd; %add gaussian noise of positional jitter
Bpos=Bpos(1:M,:); %chop off the extra neurons


%            A     B
truePairs=[mapping', (1:N)' ];
truePairs(truePairs(:,2)>M,:)=[];

%% Plot simulated position data

%Plot Acquisition A
figh=figure;
h(1)=subplot(2,2,1);
plot3(Apos(:,1),Apos(:,2),Apos(:,3),Amarker)
axis equal %set each axis to have equal spaced scaling in the plot 
axis vis3d %preserve scaling during rotations
title('Acquisition A')
xlabel('x (microns)')
ylabel('y (microns)')
zlabel('z (microns)')


%Plot Acquisition B
h(2)=subplot(2,2,2);
plot3(Bpos(:,1),Bpos(:,2),Bpos(:,3),Bmarker)
axis equal %set each axis to have equal spaced scaling in the plot 
axis vis3d %preserve scaling during rotations
title('Acquisition B')
xlabel('x (microns)')
ylabel('y (microns)')
zlabel('z (microns)')


%% Calculate Distances between all possible paris of neurons between A & B

%We are going to build up an M rows by N column matrix whereby each element
%is the euclidian distance between the neurons inexed by that row and
%column.
%
% So the element at the m'th row and n'th  column corresponds to 
% the distance between the m'th neuron from B and the n'th neuron from A

%Insead of making a bunch of boring for-loops, we are going to implement
%this in a classy and clever way.

% Here is the theory first. We note that A and B already have neurons listed
% in some sort of order. If we take the given order, it is easy to find 
% euclidean distnces between each neuron in A and its corresponding neuron
% in B. This is:
%  D= sqrt(   sum( (Apos(1:M,:)-Bpos(1:M,:)).^2 ,2 ) )
% Note we are excluding those neurons from B that have no partner in A,
% recall that N>=M.
%
% Now if we were to merely permute the order of A systematically, (by say a
% circular shift) then we
% would be able to calculate the euclidean distance between every possible
% pairing between a neuron from A and B. 
%
% We store this information in the form of a matrix.

%initialize the distance matrix
Dmat=zeros(M,N); %rows are neurons from B and columns are neurons from A

%calculate all of the distances
for k=0:N-1 %For all neurons in A
   Atemp=circshift(Apos,-k); %circullary shift the order of neurons in A
   Aind=circshift((1:N)',-k); %do the same to an array of indices to keep track of what the indices are up to
   D= sqrt(  sum( (Atemp(1:M,:)-Bpos(1:M,:)).^2 ,2)  ); % calculate the distance
   
   %get element # from row,col:  B       A                        
   indices=sub2ind(size(Dmat),[1:M]',[Aind(1:M)]); %figure out which diaganol in the distance matrix we are on
   
   Dmat(indices)=D; %record the distances on the appropriate diagonal
end

%plot distance matrix
figure;
imagesc(Dmat);
title('Distance matrix')


trueEnergy=sum(Dmat(sub2ind(size(Dmat),truePairs(:,2),truePairs(:,1))).^2);

%Check to see if the inferredPairs match the true Pairs


%% test Ben's code
Allpos=[Apos ; Bpos];

Allpos=[Allpos,[1:N,mapping(1:M)]'];
alltime=zeros(length(Apos)+length(Bpos),1);
alltime(length(Apos)+1:end)=1;
Allpos=[Allpos alltime];
param.excessive=1;
param.dims=3;
try
tracks = track(Allpos,3.5*pstd,param);
catch
    tracks = track(Allpos,2*pstd,param);
end

%% Dispaly partnered neurons


% 
% %Get ready to plot lines connected the inferred and true pairings.
% %Create 3 matrices, one for each of x,y,z where each row has the
% %coordinates for each pair of neurons
% AA=1; BB=2; XX=1; YY=2; ZZ=3;
% linesX=[Apos(inferredPairs(:,AA),XX) Bpos(inferredPairs(:,BB),XX)];
% linesY=[Apos(inferredPairs(:,AA),YY) Bpos(inferredPairs(:,BB),YY)];
% linesZ=[Apos(inferredPairs(:,AA),ZZ) Bpos(inferredPairs(:,BB),ZZ)];
% 
% linesTrueX=[Apos(truePairs(:,AA),XX) Bpos(truePairs(:,BB),XX)];
% linesTrueY=[Apos(truePairs(:,AA),YY) Bpos(truePairs(:,BB),YY)];
% linesTrueZ=[Apos(truePairs(:,AA),ZZ) Bpos(truePairs(:,BB),ZZ)];
% 
% 
% 
% figure(figh)
% h(3)=subplot(2,2,3);
% hold on;
% plot3(Apos(:,1),Apos(:,2),Apos(:,3),Amarker)
% plot3(Bpos(:,1),Bpos(:,2),Bpos(:,3),Bmarker)
% plot3(linesX',linesY',linesZ','g'); %inferred Pairs
% 
% 
% 
% axis equal %set each axis to have equal spaced scaling in the plot 
% axis vis3d %preserve scaling during rotations
% title('Inferred Pairings')
% xlabel('x (microns)')
% ylabel('y (microns)')
% zlabel('z (microns)')
% 
% 
% 
% 
% figure(figh)
% h(4)=subplot(2,2,4);
% hold on;
% plot3(Apos(:,1),Apos(:,2),Apos(:,3),Amarker)
% plot3(Bpos(:,1),Bpos(:,2),Bpos(:,3),Bmarker)
% 
% 
% plot3(linesTrueX',linesTrueY',linesTrueZ','m'); %truePairs 
% axis equal %set each axis to have equal spaced scaling in the plot 
% axis vis3d %preserve scaling during rotations
% title('True Pairings')
% xlabel('x (microns)')
% ylabel('y (microns)')
% zlabel('z (microns)')



% Make it so that when you rotate one subplot, the other rotates also
% Taken from this example: edit ([docroot '/techdoc/ref/examples/doc_linkprop']);
% Link the CameraPosition and CameraUpVector properties of each subplot axes
% hlink = linkprop(h,{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle'});
% key = 'graphics_linkprop';
% % Store link object on first subplot axes
% setappdata(h(1),key,hlink); 
% 
% 
% %display the distance  between neurons for each inferred Pair
% dist=Dmat(sub2ind(size(Dmat),inferredPairs(:,2),inferredPairs(:,1) ));
% cumulativeDist=sum(dist);
numFailsBen=M-sum(diff(tracks(:,4))==0);
disp(['Of ' num2str(M) ' possible neuron pairs,']); 
disp([num2str(M-numFailsBen) ' pairs were correctly identified by track.'])

%disp([num2str(numFails) ' pairs were incorrect.'])

disp(['true energy: ' num2str(trueEnergy)])
%disp(['final energy: ' num2str(DEnergyMin)])


%% make position matrix..try clustering?
% 
% Adis=pdist2(Apos,Apos);
% Bdis=pdist2(Bpos,Bpos);
% 
% 
% idxNames=mapping(1,:);
% Acluster=clustergram(Adis);
% Bcluster=clustergram(Bdis);
% Bidx=get(Bcluster,'RowLabels');
% Bidx=str2double(Bidx);
% Bdis=Bdis(Bidx,Bidx);
% Aidx=get(Acluster,'RowLabels');
% Aidx=str2double(Aidx);
% Adis=Adis(Aidx,Aidx);


%% sort position matricies, try to look for similarities looking within some window after sorting


Adis=pdist2(Apos,Apos);
Bdis=pdist2(Bpos,Bpos);

[~,Aidx]=sort(max(Adis));
Asort2=(Adis(Aidx,Aidx));
Asort=sort(Asort2);

[~,Bidx]=sort(max(Bdis));
Bsort2=Bdis(Bidx,Bidx);
Bsort=sort(Bsort2);

Atri=tripletTriangleAreas(Apos);
Atri2=reshape(Atri,N,[]);
Btri=tripletTriangleAreas(Bpos);
Btri2=reshape(Btri,M,[]);
tribins=linspace(0,max(Atri2(:)),200);
xbins=linspace(0,max(Asort(:)),200);
Ahist=[];
Bhist=[];
AtriHist=[];
BtriHist=[];
for iPnt=1:N
    Ahist(:,iPnt)=hist(Asort(:,iPnt),xbins);
    AtriHist(:,iPnt)=hist(Atri2(iPnt,:),tribins) ;
end
for iPnt=1:N
    Bhist(:,iPnt)=hist(Bsort(:,iPnt),xbins);

    BtriHist(:,iPnt)=hist(Btri2(iPnt,:),tribins) ;
end

ACDF=cumsum(Ahist);
BCDF=cumsum(Bhist);
ACDF=ACDF/max(ACDF(:));
BCDF=BCDF/max(BCDF(:));


AtriCDF=cumsum(AtriHist);
AtriCDF=AtriCDF/max(AtriCDF(:));
AtriCDF=AtriCDF(:,Aidx);
BtriCDF=cumsum(BtriHist);
BtriCDF=BtriCDF/max(BtriCDF(:));
BtriCDF=BtriCDF(:,Bidx);

AtriPDF=diff(AtriCDF);
APDF=diff(ACDF);
BtriPDF=diff(BtriCDF);
BPDF=diff(BCDF);

 Asort=[AtriCDF;ACDF];
 Bsort=[BtriCDF;BCDF];

 Asort=[AtriPDF; APDF];
 Bsort=[BtriPDF; BPDF];

%%
[~,unAidx]=sort(Aidx);
[~,unBidx]=sort(Bidx);
AsortTrue=Asort(:,unAidx);
AsortTrue=AsortTrue(:,mapping);
AsortTrue=AsortTrue(:,Bidx);


%%

numPoints=size(Bsort,1);
iSort=1:N;
Asorted=Asort;
%simple switching search
for iColumn=1:length(Bidx);

    tempdiff=bsxfun(@minus, Asorted,Bsort(:,iColumn));
    columnDistance=sqrt(sum(tempdiff.^2));
    [~,minDistance]=min(columnDistance);
    
    
    Asorted(:,[iColumn,minDistance])=Asorted(:,[minDistance,iColumn]);
    iSort([iColumn,minDistance])=iSort([minDistance,iColumn]);
    
end


%brute force all permutations within some window
for iColumn=1:length(Bidx);
tempIdx=max(1,iColumn-3):min(iColumn+3,length(Bidx));
    tempIdxs=permute(perms(tempIdx),[3,2,1]);
    Atemps=Asorted(:,tempIdxs);
    Atemps=reshape(Atemps,numPoints,length(tempIdx),[]);
    
    tempdiff=bsxfun(@minus, Atemps,Bsort(:,tempIdx));
    
    columnDistance=sum(sqrt(sum(tempdiff.^2)));
    [~,minDistance]=min(columnDistance);
    
    
    Asorted(:,tempIdx)=Asorted(:,tempIdxs(:,:,minDistance));
    iSort(tempIdx)=iSort(tempIdxs(:,:,minDistance));
    
end
%%

%so something markovy?
Aenergy=sqrt(sum((Asorted-Bsort).^2));
Aenergy=.5*(Asorted.*log(Asorted./(.5*(Asorted+Bsort)))+Bsort.*log(Bsort./(.5*(Asorted+Bsort))));
Aenergy=nansum(Aenergy);
Amean=mean(Aenergy);
pointList=1:length(Bidx);
tic
for iColumn=1:length(Bidx);
    for iIteration=1:150
    tempIdx=max(1,iColumn-25):min(iColumn+25,length(Bidx));

    tempIdx=tempIdx(tempIdx~=iColumn);
    otherE=Aenergy(tempIdx);
    tempIdx=randsample(tempIdx,4,1,exp(otherE/mean(otherE)));
tempIdx=[iColumn,unique(tempIdx)];
%tempIdx=tempIdx(Aenergy(tempIdx)>(2*Amean));
    tempIdxs=permute(perms(tempIdx),[3,2,1]);
    Atemps=Asorted(:,tempIdxs);
    Atemps=reshape(Atemps,numPoints,length(tempIdx),[]);
    
  %  tempdiff=bsxfun(@minus, Atemps,Bsort(:,tempIdx));
    
        avgSort=bsxfun(@plus, Atemps,Bsort(:,tempIdx))/2;
        
        JS1=.5*Atemps.*log(Atemps./avgSort);
        JS2=bsxfun(@rdivide, Bsort(:,tempIdx),avgSort);
        JS2=.5*bsxfun(@times, Bsort(:,tempIdx), log(JS2));
        tempdiff=(JS1+JS2)/2;

%     for iPerm=1:size(Atemps,3)
%         columnDistance(iPerm)=corr2(Atemps(:,:,iPerm),Bsort(:,tempIdx));
%     end
    
    
   columnDistance=sum(sqrt(sum(tempdiff.^2)));
    [~,minDistance]=min(columnDistance);
    
    
    Asorted(:,tempIdx)=Asorted(:,tempIdxs(:,:,minDistance));
    iSort(tempIdx)=iSort(tempIdxs(:,:,minDistance));
    Aenergy=sqrt(sum((Asorted-Bsort).^2));
    Aenergy=.5*(Asorted.*log(Asorted./(.5*(Asorted+Bsort)))+Bsort.*log(Bsort./(.5*(Asorted+Bsort))));
Aenergy=nansum(Aenergy);
Amean=mean(Aenergy);

    end
end
toc

numberCorrect=sum(unmapping(Aidx(iSort))'==Bidx');
disp([num2str(sum(unmapping(Aidx(iSort))'==Bidx')) ' pairs were correctly identified by distance metrics.'])
AenergyTrue=.5*(AsortTrue.*log(AsortTrue./(.5*(AsortTrue+Bsort)))+Bsort.*log(Bsort./(.5*(AsortTrue+Bsort))));
AenergyTrue=nansum(Aenergy);
AenergyTrue=(sqrt(sum((AsortTrue-Bsort).^2)));
    Aenergy=sqrt(sum((Asorted-Bsort).^2));

dmatEnergy=Apos((Aidx(iSort))',:)-Bpos(Bidx',:);
dmatEnergy=sum((sum(dmatEnergy.^2,2)));

disp(['final dmat energy: ' num2str(dmatEnergy)])

%%
Asort2=Adis(Aidx(iSort),Aidx(iSort));

imagesc(~bsxfun(@or, unmapping(Aidx(iSort))'==Bidx',unmapping(Aidx(iSort))==Bidx).*Asort2)

Aenergy=sqrt(sum((Asort2-Bsort2).^2));
Amean=mean(Aenergy);
pointList=1:length(Bidx);
tic
for iColumn=1:length(Bidx);
    for iIteration=1:150
    tempIdx=max(1,iColumn-25):min(iColumn+25,length(Bidx));

    tempIdx=tempIdx(tempIdx~=iColumn);
    otherE=Aenergy(tempIdx);
    tempIdx=randsample(tempIdx,4,1,exp(otherE/mean(otherE)));
tempIdx=[iColumn,unique(tempIdx)];
%tempIdx=tempIdx(Aenergy(tempIdx)>(2*Amean));
tempIdxs2=perms(tempIdx);
    tempIdxs=permute(tempIdxs2,[3,2,1]);
    tempIdxs2=permute(tempIdxs2,[2,3,1]);
    
    Atemps=Asort2(:,tempIdxs);
    AAtemp=Asort2(tempIdxs,:);
    Atemps=reshape(Atemps,numPoints,length(tempIdx),[]);
    AAtemp=permute(reshape(AAtemp,length(tempIdx),[],numPoints),[1,3,2]);

        % fix the indexing
        
        for i=1:size(Atemps,3)
    Atemps(tempIdx,:,i)=AAtemp(:,tempIdxs2(:,:,i),i);
        end
        
    tempdiff=bsxfun(@minus, Atemps,Bsort(:,tempIdx));
    
    columnDistance=sum(sqrt(sum(tempdiff.^2)));
    [~,minDistance]=min(columnDistance);
    
    
    Asort2(:,tempIdx)=Asort2(:,tempIdxs(:,:,minDistance));
        Asort2(tempIdx,:)=Asort2(tempIdxs(:,:,minDistance),:);

    iSort(tempIdx)=iSort(tempIdxs(:,:,minDistance));
    Aenergy=sqrt(sum((Asort2-Bsort2).^2));
Amean=mean(Aenergy);
    end
end
toc
