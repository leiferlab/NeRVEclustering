
totalPoints=20;
fiducialFiles=uipickfiles('filterspec','O:\');
%%

    data=load(fiducialFiles{1});
    fpoints=data.fiducialPoints;
    maxL=max(cell2mat(cellfun(@(x) length(x), fpoints,'uniform',0)));
    dataAll=cell(totalPoints,5);
    dataAll=repmat({dataAll}, length(fpoints),1);
    nTimes=length(fpoints);
  %%
    start=1;
   startVec=[];
for iFiles=1:length(fiducialFiles);
    data=load(fiducialFiles{iFiles});
    fpoints=data.fiducialPoints;
    L=cell2mat(cellfun(@(x) size(cell2mat(x(:,1:4)),1), fpoints,'uniform',0));
    maxL=max(L);
if maxL>0
    for iTime=1:nTimes
        dataAll{iTime}(start:start+L(iTime)-1,:)=fpoints{iTime}(1:L(iTime),:);
    end
    
end
startVec=[startVec,start-1];
    start=start+maxL;
   [~, tempName]=fileparts(fiducialFiles{iFiles});
   FileName{iFiles}=tempName;
end
%%
dataFolder=fileparts(fiducialFiles{1});
dataFolder=fileparts(dataFolder);
fiducialPoints=dataAll;
timeOffset=-6;

save([dataFolder filesep 'tempFiducials 20150327'],'fiducialPoints','timeOffset');

%%
dataFolder='O:\20150118\BrainScanner20150118_184857 - 4+min';
load([dataFolder filesep 'tempFiducials20140129'],'timeOffset','fiducialPoints');
dataAll=fiducialPoints;

%%
xAll=[];tAll=[];xMatAll=[];
hasPoints=cellfun(@(x) ~isempty(x{1}), dataAll,'uniformoutput',0);
hasPoints=find(cell2mat(hasPoints));
nTimes=length(hasPoints);


for ii=1:nTimes
    i=hasPoints(ii);
%x=(cell2mat(cellfun(@(x) x{:,1}, fiducialPoints(i),'uniformoutput',0)));
x=(dataAll{i});



emptyX=cellfun(@(x) ~isempty(x),x(:,1));
plotIdx=find(emptyX);
if ii==1
X0=zeros(max(plotIdx));
end

x=cell2mat(x);
x=x(:,1:2);
X2=(pdist(x));
X2=squareform(X2);

max(plotIdx)
Xtemp=X0;
Xtemp(plotIdx,plotIdx)=X2;
X2=squareform(Xtemp);

    xAll=cat(1,xAll,X2);
    xMatAll=cat(3,xMatAll,squareform(X2));
    tAll=cat(1,tAll,i);



end
%%
xAll(xAll==0)=nan;
xmeans=nanmean(xAll,1);
xSTDs=nanstd(xAll,[],1);
zAll=bsxfun(@minus, xAll,xmeans);
zAll=bsxfun(@rdivide, zAll,xSTDs);

[coeff,score,xlatent,tsquared,explained,mu] = pca(zAll,'row','pairwise');
coeff2=bsxfun(@times, coeff,xSTDs');
coeff2=bsxfun(@plus, coeff2,xmeans');


%%
connAll=false(totalPoints,totalPoints,length(dataAll));
for iTime=1:length(dataAll);
t=dataAll{iTime};
emptyT=cellfun(@(x) ~isempty(x),t(:,1));
plotIdx=find(emptyT);
t2=cell2mat(t);




if ~isempty(t2)
tsize=size(t2,1);
TRI=delaunayTriangulation(t2(:,1),t2(:,2),t2(:,4));
conn=TRI.ConnectivityList;

edges1=edges(TRI);
edgeIdx=sub2ind([tsize tsize], edges1(:,1),edges1(:,2));
connMatrix=false(tsize);
connMatrix(edgeIdx)=true;
connMatrix=connMatrix | connMatrix';
connFull=false(length(t));
connFull(plotIdx,plotIdx)=connMatrix;
connAll(:,:,iTime)=connFull;
end
end
%%
connAll=double(connAll);
zConnAll=zscore(connAll,0,3);
