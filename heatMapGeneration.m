
function heatMapGeneration(dataFolder,runTime)


% take a datafolder, loading the wormFiducialIntensities to produce
% heatmaps and ordering for Red, green and ratios
if nargin==0
    dataFolder=uipickfiles();
    dataFolder=dataFolder{1};
end
%%
fiducialIFiles=dir([dataFolder filesep 'wormFiducial*']);
fiducialIFiles={fiducialIFiles.name};
pointFiles=cellfun(@(x) strfind(x,'oint'),fiducialIFiles,'Uniform',0);
pointFiles=cellfun(@(x) ~isempty(x), pointFiles);
pointFiles=find(pointFiles);
%%
if ~isempty(pointFiles)
load([dataFolder filesep fiducialIFiles{pointFiles(1)}])
RvalsAll=RvalAll;
GvalsAll=GvalAll;

else
load([dataFolder filesep fiducialIFiles{1}])
ballSearch=7;
cropSize=size(Rcrop);
nTimes=cropSize(4);
vols=cropSize(5);
RvalsAll=zeros(vols,nTimes);
GvalsAll=RvalsAll;

center=ceil(cropSize/2);
center=center(1:3);
centerIdx=sub2ind(cropSize(1:3),center(1),center(2),center(3));
[ballx,bally]=meshgrid(1:cropSize(1),1:cropSize(2));
ballx=ballx-center(1);
bally=bally-center(2);
ballr=sqrt(ballx.^2+bally.^2);
ballr=repmat(ballr<ballSearch,1,1,cropSize(3));
progressbar(0,0);
for iVol=1:vols
    for i=1:nTimes
      %  iTime=hasPoints(i);
      iTime=i;
        progressbar(iVol/vols,iTime/nTimes);
        subVolumeR=(Rcrop(:,:,:,iTime,iVol))-90;
     %           subVolumeR=pedistalSubtract(Rcrop(:,:,:,iTime,iVol));

        subVolumeG=(Gcrop(:,:,:,iTime,iVol))-90;
        subVolumeG((subVolumeG)<=0)=nan;
        subVolumeR(subVolumeR<0)=nan;

        subBW=(subVolumeR)>(.5*subVolumeR(center(1),center(2),center(3)));
        subBW=subBW.*ballr;
        cc=bwconncomp(subBW,6);
        centerPixPresent=cellfun(@(x) any(x==centerIdx),cc.PixelIdxList);
        if any(centerPixPresent)
            cc.PixelIdxList=cc.PixelIdxList(centerPixPresent);
            areaIdx=1;
        else
        
        nPixRegions=cell2mat(cellfun(@(x) length(x), cc.PixelIdxList,'uniformOutput',0));
        [~,areaIdx]=max(nPixRegions);
        end
        Rvals=subVolumeR(cell2mat(cc.PixelIdxList(areaIdx)'));
        Gvals=subVolumeG(cell2mat(cc.PixelIdxList(areaIdx)'));
        
        RvalsAll(iVol,i)=nanmean(Rvals);
        GvalsAll(iVol,i)=nanmean(Gvals);
    end
end

%%
       RvalsAll(RvalsAll==0|RvalsAll>(nanstd(RvalsAll(:))*10))=nan;

GvalsAll(GvalsAll==0|GvalsAll>(nanstd(GvalsAll(:))*10))=nan;
end
% rejects=(mean(isnan(RvalsAll),2)>.1);
% rejectCorrection=cumsum(rejects);
% rejectCorrection(rejects)=[];
% rejectCorrection=rejectCorrection+(1:length(rejectCorrection))';
% GvalsAll(rejects,:)=[];
% 
% RvalsAll(rejects,:)=[];
%%

%A=(G2smooth./R2smooth)';
%correct photobleaching red
photoBleachingR=zeros(size(RvalsAll));
photoBleachingG=zeros(size(RvalsAll));

%Fexponent =@(x,xdata) x(1)*exp(x(2)*xdata)+x(3);
Fexponent=fittype('a*exp(b*x)+c','dependent',{'y'},'independent',...
    {'x'},'coefficients',{'a', 'b', 'c'});
fitOptions=fitoptions(Fexponent);
fitOptions.Lower=[0,-.2,0];
fitOptions.Upper=[1000,0,10000];

%progressbar(0)
for i=1:size(RvalsAll,1)
    
    try
        %%
     
%    progressbar(i/size(RvalsAll,1));
    xVals=(1:size(RvalsAll,2))';
    
    present=(~isnan(RvalsAll(i,:)+GvalsAll(i,:))') ;
    xVals=xVals(present);
    rVals=RvalsAll(i,:)';
    minWindow=150;
    gVals=GvalsAll(i,:)';
        gVals=gVals(present);
    rVals=rVals(present);
    gVals=ordfilt2(gVals,30,true(minWindow,1));
    rVals=ordfilt2(rVals,30,true(minWindow,1));
    gVals(1:(minWindow/2))=gVals((minWindow/2));
    rVals(1:(minWindow/2))=rVals((minWindow/2));
    gVals(end-(minWindow/2):end)=gVals(end-(minWindow/2));
    rVals(end-(minWindow/2):end)=rVals(end-(minWindow/2));
       fitOptions.StartPoint=[1,-.1,min(rVals)];


f=fit(xVals,rVals,Fexponent,fitOptions);
   fitOptions.StartPoint=[1,-.1,min(gVals)];

g=fit(xVals,gVals,Fexponent,fitOptions);
if f(1)>(max(RvalsAll(i,:))+100)
    f=fit(xVals,rVals,'poly1');
if f.p1>0
    f.p1=0;
end
    
end
if g(1)>(max(GvalsAll(i,:))+100)
    g=fit(xVals,gVals,'poly1');
if g.p1>0
    g.p1=0;
end

  %  g.a=0;
end
%[f.b g.b]
if 1
subplot(2,1,1);
plot(GvalsAll(i,:))
hold on
plot(g)
%imagesc([],[0 1000],ethoTrack')
%alpha(.2);
ylim([0 max(GvalsAll(i,:))]);
hold off
subplot(2,1,2);
plot(RvalsAll(i,:))
hold on
plot(f)
%imagesc([],[0 1000],ethoTrack')
%alpha(.2);
ylim([0 max(RvalsAll(i,:))]);
hold off
drawnow
pause(.1)
end
photoBleachingR(i,:)=f((1:size(RvalsAll,2)))-f(size(RvalsAll,2));

photoBleachingG(i,:)=g((1:size(RvalsAll,2)))-g(size(RvalsAll,2));
    catch me
        me
    end
    

end
%%
% photoBleachingR=bsxfun(@minus,photoBleachingR,photoBleachingR(:,end));
% photoBleachingG=bsxfun(@minus,photoBleachingG,photoBleachingG(:,end));
Rvalstemp=RvalsAll-photoBleachingR ;

Rvalstemp(bsxfun(@le,Rvalstemp,quantile(Rvalstemp,.1,2)))=nan;
Rvalstemp=colNanFill(Rvalstemp')';

Gvalstemp=GvalsAll-photoBleachingG ;

Gvalstemp(bsxfun(@le,Gvalstemp,quantile(Gvalstemp,.1,2)))=nan;
Gvalstemp=colNanFill(Gvalstemp')';


%%
%f = lsqcurvefit(Fexponent,[.001,-.001,min(photoBleaching)/2],hasPointsTime,photoBleaching')
A=Rvalstemp';

A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=medfilt2(A2,[3,1]);
A2=smooth2a(A2,2,0);
A(~isnan(A))=A2(~isnan(A));
A(A<-1)=-nan;
R2=A';


A=(Gvalstemp)';
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=medfilt2(A2,[3,1]);
A2=smooth2a(A2,2,0);
A(~isnan(A))=A2(~isnan(A));
A(A<-1)=-nan;
G2=A';
%A=medfilt2(A,[3,1]);
%A=R2smooth';
% A(R2smooth<.1)=nan;
% A=[A score(1:size(A,1),1:3)];
% A=G2smooth';

A=(Gvalstemp)'./Rvalstemp';
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
%A(~isnan(A))=A2(~isnan(A));
%A(A<-1)=-nan;
A2=colNanFill(A);
A2=medfilt2(A2,[3,1]);
A2=smooth2a(A2,2,0);

Ratio2=A2';

%A=smooth2a(A,0,1);
A(isnan(A))=0;
acorr=corr(A);
atemp=nancov(A)./sqrt(nanvar(A)'*nanvar(A));
acorr(isnan(acorr))=atemp(isnan(acorr));
acorr(isnan(acorr))=0;
% dcorr=1-acorr;
% dcorr(eye(size(dcorr))>0)=0;
% cg=linkage(squareform(dcorr));
% cg2=cluster(cg,.001);
% cg2=cg2+(1:length(cg2))'/1000;
% [~,ia]=sort(cg2);
% [~,ib]=sort(ia);

cg = clustergram(acorr);
cgIdx=str2double(get(cg,'RowLabels'));
[~,cgIdxRev]=sort(cgIdx);
%cgIdx=cgIdx(cgIdx<(max(cgIdx)-2));







%%
save([dataFolder filesep 'heatData' runTime],'G2','R2','Ratio2','acorr','cgIdx','cgIdxRev');

