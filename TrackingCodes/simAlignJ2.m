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

N=125; %acquisition A
M=120; %acquisition B
assert(N>=M); %by convention, the code assumes A has more neurons then B

%set dimensions of  simulated 3D volume in microns
xdiam=75;
ydiam=150;
zdiam=75;

%set standard deviation of jitter in neural position between acquisition A & B
pstd=6;


%Set the markers to use when plotting A & B
Bmarker='rs';
Amarker='bo';

%% Simulate position data 

Apos=[rand(N,1).*xdiam, rand(N,1).*ydiam, rand(N,1).*zdiam];

mapping=randperm(N);  %Elements are neuron ID's from A  and their index is the neuorn ID from B
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

%% Partner up all of the neurons, one in A and one in B

% The algorithm we will use is the following:
% First find the pairing that minimizes the distance
% Then eliminate both neurons from future consideration
% Repeat.
%
% Note we will have to be clever in our bookkeeping.
Dweight=exp(-Dmat.^2/2/pstd);
%Dweight=1./(1+Dmat.^2/2/pstd^3);
% imagesc(Dweight);
% hold on
% p=[];
% scatter(truePairs(:,1),truePairs(:,2),'xr');

%Dweight=Dweight(:);
[colIDs,rowIDs]=meshgrid(1:N,1:M);
rowIDs=rowIDs(:);
colIDs=colIDs(:);
data=1:numel(M);
DEnergyMin=Inf;
pairsBest=zeros(M,2);
linidxmap=zeros(M,500);
for iterations=1:500;

Dtemp=Dweight;
pairs=zeros(M,2);
linidx=zeros(M,1);
rowSearch=randperm(M);

for k=1:M % there will be M pairs (recall N>=M)
    %draw random row
    %draw term from distribution
    [~,idx] = datasample(Dweight(:),1,'Weights',Dtemp(:)) ;
    [m,n]=ind2sub(size(Dmat),idx);

    pairs(k,:)=[n,m];
    %Exclude all elements in row m and all elemnts in column n from future
    %consideration because we have already paired those neurons
    Dtemp(rowIDs==m)=0;
    Dtemp(colIDs==n)=0;
end
% linidx=sort(linidx);
% if any(all(bsxfun(@eq,linidxmap,linidx),1))
% display('old path');
% end
% linidxmap(:,iterations)=linidx;

DEnergy=sum(Dmat(sub2ind(size(Dmat),pairs(:,2),pairs(:,1))).^2);
if DEnergy<DEnergyMin
    DEnergyMin=DEnergy;
    pairsBest=pairs;
 
end
%    if ~isempty(p) && all(ishandle(p)), delete(p); end
% 
%     p=scatter(pairs(:,1),pairs(:,2),'og');
% pause(.01);
%     
end

pairs=pairsBest;
inferredPairs=sortrows(pairs,2);
truePairs;
trueEnergy=sum(Dmat(sub2ind(size(Dmat),truePairs(:,2),truePairs(:,1))).^2);

%Check to see if the inferredPairs match the true Pairs

difference=inferredPairs-truePairs;


%% Dispaly partnered neurons



%Get ready to plot lines connected the inferred and true pairings.
%Create 3 matrices, one for each of x,y,z where each row has the
%coordinates for each pair of neurons
AA=1; BB=2; XX=1; YY=2; ZZ=3;
linesX=[Apos(inferredPairs(:,AA),XX) Bpos(inferredPairs(:,BB),XX)];
linesY=[Apos(inferredPairs(:,AA),YY) Bpos(inferredPairs(:,BB),YY)];
linesZ=[Apos(inferredPairs(:,AA),ZZ) Bpos(inferredPairs(:,BB),ZZ)];

linesTrueX=[Apos(truePairs(:,AA),XX) Bpos(truePairs(:,BB),XX)];
linesTrueY=[Apos(truePairs(:,AA),YY) Bpos(truePairs(:,BB),YY)];
linesTrueZ=[Apos(truePairs(:,AA),ZZ) Bpos(truePairs(:,BB),ZZ)];



figure(figh)
h(3)=subplot(2,2,3);
hold on;
plot3(Apos(:,1),Apos(:,2),Apos(:,3),Amarker)
plot3(Bpos(:,1),Bpos(:,2),Bpos(:,3),Bmarker)
plot3(linesX',linesY',linesZ','g'); %inferred Pairs



axis equal %set each axis to have equal spaced scaling in the plot 
axis vis3d %preserve scaling during rotations
title('Inferred Pairings')
xlabel('x (microns)')
ylabel('y (microns)')
zlabel('z (microns)')




figure(figh)
h(4)=subplot(2,2,4);
hold on;
plot3(Apos(:,1),Apos(:,2),Apos(:,3),Amarker)
plot3(Bpos(:,1),Bpos(:,2),Bpos(:,3),Bmarker)


plot3(linesTrueX',linesTrueY',linesTrueZ','m'); %truePairs 
axis equal %set each axis to have equal spaced scaling in the plot 
axis vis3d %preserve scaling during rotations
title('True Pairings')
xlabel('x (microns)')
ylabel('y (microns)')
zlabel('z (microns)')



% Make it so that when you rotate one subplot, the other rotates also
% Taken from this example: edit ([docroot '/techdoc/ref/examples/doc_linkprop']);
% Link the CameraPosition and CameraUpVector properties of each subplot axes
hlink = linkprop(h,{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle'});
key = 'graphics_linkprop';
% Store link object on first subplot axes
setappdata(h(1),key,hlink); 


%display the distance  between neurons for each inferred Pair
dist=Dmat(sub2ind(size(Dmat),inferredPairs(:,2),inferredPairs(:,1) ));
cumulativeDist=sum(dist);

numFails=sum(difference(:,1)~=0);
disp(['Of ' num2str(M) ' possible neuron pairs,']); 
disp([num2str(M-numFails) ' pairs were correctly identified.'])
disp([num2str(numFails) ' pairs were incorrect.'])
disp(['true energy: ' num2str(trueEnergy)])
disp(['final energy: ' num2str(DEnergyMin)])