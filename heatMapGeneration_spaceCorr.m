
function heatMapGeneration_spaceCorr(dataFolder,fiducialPoints,alignments,runTime)


% take a datafolder, loading the wormFiducialIntensities to produce
% heatmaps and ordering for Red, green and ratios
if nargin==0
    dataFolder=uipickfiles();
    dataFolder=dataFolder{1};
end

%%
fiducialIFiles=dir([dataFolder filesep 'WormFiducial*']);
fiducialIFiles={fiducialIFiles.name};
pointFiles=cellfun(@(x) strfind(x,'oint'),fiducialIFiles,'Uniform',0);
pointFiles=cellfun(@(x) ~isempty(x), pointFiles);
pointFiles=find(pointFiles);
%%
if ~isempty(pointFiles)
    if length(pointFiles)>1
        display('Select point files');
        pointFiles=uipickfiles('filterspec',dataFolder);
        pointFiles=pointFiles{1};
        load(pointFiles)
    else
    load([dataFolder filesep fiducialIFiles{pointFiles(1)}])

    end
RvalsAll=RvalAll;
GvalsAll=GvalAll;
RvalsAll=RvalsAll-90;
GvalsAll=GvalAll-90;
RvalsAll(RvalsAll<25)=nan;
GvalsAll(isnan(RvalsAll))=nan;
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




[FX,FY,~]=fiducials2mat(fiducialPoints,size(GvalsAll,1));
%% lower correction
% FR=(FX-nanmean(FX(:))).^2+(FY-nanmean(FY(:))).^2;
% FR=sqrt(FR);
load('IlluminationProfile.mat')
rect1=alignments.S2AHiRes.rect1;
rect2=alignments.S2AHiRes.rect2;


upper=centroids(:,2)>650;
lower=centroids(:,2)<550;
cen_upper=centroids(upper,:);
cen_lower=centroids(lower,:);
I_lower=intensities(lower);
I_upper=intensities(upper);

cen_upper=bsxfun(@minus,cen_upper,rect2(1:2)-1);
cen_lower=bsxfun(@minus,cen_lower,rect1(1:2)-1);

cen_upper=transformPointsForward(alignments.S2AHiRes.t_concord,cen_upper);


fupper=fit(cen_upper,I_upper,'poly33');
flower=fit(cen_lower,I_lower,'poly33');

FX(FX<5)=nan;
FY(FY<5)=nan;
FX=colNanFill(FX')';
FY=colNanFill(FY')';
Rcorr=flower(FY-rect1(1)+1,FX-rect1(2)+1);
Rcorr=Rcorr/nanmean(Rcorr(:));
Gcorr=fupper(FY-rect1(1)+1,FX-rect1(2)+1);
Gcorr=Gcorr/nanmean(Gcorr(:));
Rcorr(isnan(Rcorr))=1;
Gcorr(isnan(Gcorr))=1;

% rejects=(mean(isnan(RvalsAll),2)>.1);
% rejectCorrection=cumsum(rejects);
% rejectCorrection(rejects)=[];
% rejectCorrection=rejectCorrection+(1:length(rejectCorrection))';
 GvalsAll2=GvalsAll./Gcorr;
% 
 RvalsAll2=RvalsAll./Rcorr;
%%

%A=(G2smooth./R2smooth)';
%correct photobleaching red
photoBleachingR=zeros(size(RvalsAll2));
photoBleachingG=zeros(size(RvalsAll2));

%Fexponent =@(x,xdata) x(1)*exp(x(2)*xdata)+x(3);
Fexponent=fittype('a*exp(b*x)+c','dependent',{'y'},'independent',...
    {'x'},'coefficients',{'a', 'b', 'c'});

fitOptions=fitoptions(Fexponent);
fitOptions.Lower=[0,-.2,0];
fitOptions.Upper=[1000,0,10000];

%progressbar(0)
for i=1:size(RvalsAll2,1)
    
    try
        %%
     
%    progressbar(i/size(RvalsAll2,1));
    xVals=(1:size(RvalsAll2,2))';
    
    present=(~isnan(RvalsAll2(i,:)+GvalsAll2(i,:))') ;
    xVals=xVals(present);
    rVals=RvalsAll2(i,:)';
    minWindow=150;
    gVals=GvalsAll2(i,:)';
        gVals=gVals(present);
    rVals=rVals(present);
    gVals=ordfilt2(gVals,30,true(minWindow,1));
    rVals=ordfilt2(rVals,30,true(minWindow,1));
   
       fitOptions.StartPoint=[range(rVals(rVals~=0)),-.0006,min(rVals(rVals~=0))];
fitOptions.Weights=zeros(size(rVals));
fitOptions.Weights(minWindow:end-minWindow)=1;
[f,fout]=fit(xVals,rVals,Fexponent,fitOptions);
if fout.rsquare<.9
    logVals=log(rVals);
    logVals=logVals(rVals~=0);
    logXvals=xVals(rVals~=0);
    expFit=polyfit(logXvals,logVals,1);
    [f.b expFit(1) 1-fout.rsquare]
    
    f.a=exp(expFit(2));
    f.b=expFit(1);
    
end

   fitOptions.StartPoint=[range(gVals),-.001,min(gVals)];
fitOptions.Weights=zeros(size(gVals));
fitOptions.Weights(minWindow:end-minWindow)=1;


[~,maxPos]=max(gVals(1:300));
fitOptions.Weights(1:maxPos)=0;
[g,gout]=fit(xVals,gVals,Fexponent,fitOptions);
if f(1)>(max(RvalsAll2(i,:))+100)
    f=fit(xVals,rVals,'poly1');
if f.p1>0
    f.p1=0;
end
    
end
if g(1)>(max(GvalsAll2(i,:))+1000)
    g=fit(xVals,gVals,'poly1');
if g.p1>0
    g.p1=0;
end

  %  g.a=0;
end
%[f.b g.b]
if 0
subplot(2,1,1);
plot(GvalsAll2(i,:))
hold on
plot(g)
ylim([0 g(0)+100])

%imagesc([],[0 1000],ethoTrack')
%alpha(.2);
hold off
subplot(2,1,2);
plot(RvalsAll2(i,:))
hold on

plot(f)
%imagesc([],[0 1000],ethoTrack')
%alpha(.2);
ylim([0 f(0)+100]);
hold off
drawnow
pause(.1)
end
photoBleachingR(i,:)=f((1:size(RvalsAll2,2)))-f(size(RvalsAll2,2));

photoBleachingG(i,:)=g((1:size(RvalsAll2,2)))-g(size(RvalsAll2,2));
    catch me
        me
    end
    

end
%%
% photoBleachingR=bsxfun(@minus,photoBleachingR,photoBleachingR(:,end));
% photoBleachingG=bsxfun(@minus,photoBleachingG,photoBleachingG(:,end));
Rvalstemp=RvalsAll2-photoBleachingR ;
RvalstempZ=bsxfun(@minus,Rvalstemp,nanmean(Rvalstemp,2));
RvalstempZ=bsxfun(@rdivide,RvalstempZ,nanstd(RvalstempZ,[],2));


Rvalstemp(RvalstempZ<-2|RvalstempZ>6|Rvalstemp<40)=nan;
%Rvalstemp=colNanFill(Rvalstemp')';

Gvalstemp=GvalsAll2-photoBleachingG ;
GvalstempZ=bsxfun(@minus,Gvalstemp,nanmean(Gvalstemp,2));
GvalstempZ=bsxfun(@rdivide,GvalstempZ,nanstd(GvalstempZ,[],2));

Gvalstemp(GvalstempZ>10|Gvalstemp<0)=nan;
%Gvalstemp=colNanFill(Gvalstemp')';


%%
%f = lsqcurvefit(Fexponent,[.001,-.001,min(photoBleaching)/2],hasPointsTime,photoBleaching')

%%
A=Rvalstemp';
Asmooth=smooth2a(A,50,0);
Asmooth=colNanFill(Asmooth);

A0=quantile(Asmooth,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
%A2=medfilt2(A2,[3,1]);
A2=imfilter(A2, gausswin(5,1)/sum( gausswin(5,1)));
A2(A2<-1)=-nan;
R2=A2';


A=(Gvalstemp)';
Asmooth=smooth2a(A,50,0);
Asmooth=colNanFill(Asmooth);

A0=quantile(Asmooth,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=imfilter(A2, gausswin(5,1)/sum( gausswin(5,1)));
A2(A2<-1)=-nan;
G2=A2';

%chop out flashes

rmean=nanmean(R2(:));
rstd=nanstd(R2(:));
nanmapr=R2>4|isnan(R2);%(rmean+3*rstd);
gmean=nanmean(G2(:));
gstd=nanstd(G2(:));
nanmapg=G2>4|isnan(G2);%(gmean+3*gstd);
G2(nanmapg)=nan;
Gvalstemp(nanmapg)=nan;
Rvalstemp(nanmapr)=nan;


%A=medfilt2(A,[3,1]);
%A=R2smooth';
% A(R2smooth<.1)=nan;
% A=[A score(1:size(A,1),1:3)];
% A=G2smooth';

A=(Gvalstemp)'./Rvalstemp';
Asmooth=smooth2a(A,20,0);
Asmooth=colNanFill(Asmooth);
A0=quantile(Asmooth,.2,1);

A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);

%A(~isnan(A))=A2(~isnan(A));
%A(A<-1)=-nan;
A2=colNanFill(A);
gfilt=@(x,h) imfilter(x, gausswin(h,1)/sum( gausswin(h,1)));
A2=imfilter(A2, gausswin(5,1)/sum( gausswin(5,1)));


Gsmooth=colNanFill(Gvalstemp');
Rsmooth=colNanFill(Rvalstemp');
Gsmooth=gfilt(Gsmooth,5);
Rsmooth=gfilt(Rsmooth,5);
A=Gsmooth./Rsmooth;
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);


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
%%
cg = clustergram(acorr);
cgIdx=str2double(get(cg,'RowLabels'));
[~,cgIdxRev]=sort(cgIdx);
%cgIdx=cgIdx(cgIdx<(max(cgIdx)-2));
close all hidden




rRaw=RvalsAll;
gRaw=GvalsAll;
rPhotoCorr=Rvalstemp;
gPhotoCorr=Gvalstemp;

%%
save([dataFolder filesep 'heatData0' runTime],'G2','R2','gRaw','rRaw',...
    'rPhotoCorr','gPhotoCorr','Ratio2','acorr','cgIdx','cgIdxRev','runTime');

