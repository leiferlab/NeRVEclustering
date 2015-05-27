fmax=3;%round(peakFreq(sin(behaviorThetaTrack)));
gDistr=load('Y:\CommunalCode\3dbrain\GFPcontrolSignalDistr.mat');
gDistr=gDistr.gDistr;
%%
bootN=10000;
deltabootRB=zeros(size(G2,1),bootN);
    deltabootGB=deltabootRB;
    deltabootRF=deltabootRB;
    deltabootGF=deltabootRB;
    deltabootRT=deltabootRB;
    deltabootGT=deltabootRB;
    deltabootGP=deltabootRB;
GbootCorrF=deltabootRB;
GbootCorrB=deltabootRB;
GbootCorrP=deltabootRB;
GbootCorrT=deltabootRB;

%%
progressbar(0)

for i=1:size(G2,1);
    progressbar(i/size(G2,1));
    gtemp=G2(i,:);
    rtemp=R2(i,:);
    
    behaviorTemp=ethoTrack';
    
   % behaviorTemp(250:end)=nan;
    anans=isnan(gtemp+rtemp)';
 %       behaviorTemp=cos(behaviorThetaTrack);

    gtemp=(gtemp(~anans));
    rtemp=(rtemp(~anans));
    behaviorTemp=behaviorTemp(~anans);
    gtemp=detrend(gtemp);
    behaviorF=(behaviorTemp==1);
    behaviorP=(behaviorTemp==-1|behaviorTemp==4);
    behaviorB=(behaviorTemp==-1);   
        behaviorT=(behaviorTemp==4);   

        bnans=isnan(ethoTrack);
        bnans=bnans(:);
%     zTemp=(vTrack(~(bnans|anans)));
%     gtemp2=gtemp(~(bnans|anans));
%     rtemp2=rtemp(~(bnans|anans));
%     
   %gzcorr(i)=corr2JN(gtemp2,zTemp');

  if 0
      
    thetaTemp=behaviorThetaTrack(~(bnans|anans));
    gtemp2=gtemp(~(bnans|anans));
    rtemp2=rtemp(~(bnans|anans));

      
   
    [phi(i), thetaCorr(i)]=fminsearch(@(x) -(corr2JN(rtemp2,sin(thetaTemp+x))),0);
   sineTemp=sin(thetaTemp+phi(i));
    thetaProj(i)=dot(sineTemp,rtemp2)/(sineTemp*sineTemp');
    
    R3(i,:)=R2(i,:)-sin(behaviorThetaTrack+phi(i))*thetaProj(i);
  end  
  
  
    RCorrB(i)=corr2JN(rtemp,behaviorB);
    GcorrB(i)=corr2JN(gtemp,behaviorB);
    RCorrF(i)=corr2JN(rtemp,behaviorF);
    GcorrF(i)=corr2JN(gtemp,behaviorF);
    RCorrP(i)=corr2JN(rtemp,behaviorP);
    GcorrP(i)=corr2JN(gtemp,behaviorP);
    RCorrT(i)=corr2JN(rtemp,behaviorT);
    GcorrT(i)=corr2JN(gtemp,behaviorT);
    deltaRB(i)=nanmean(rtemp(behaviorB))-nanmean(rtemp(~behaviorB));
    deltaGB(i)=nanmean(gtemp(behaviorB))-nanmean(gtemp(~behaviorB));
    deltaRF(i)=nanmean(rtemp(behaviorF))-nanmean(rtemp(~behaviorF));
    deltaGF(i)=nanmean(gtemp(behaviorF))-nanmean(gtemp(~behaviorF));
    deltaRT(i)=nanmean(rtemp(behaviorT))-nanmean(rtemp(~behaviorT));
    deltaGT(i)=nanmean(gtemp(behaviorT))-nanmean(gtemp(~behaviorT));
    deltaRP(i)=nanmean(rtemp(behaviorP))-nanmean(rtemp(~behaviorP));
    deltaGP(i)=nanmean(gtemp(behaviorP))-nanmean(gtemp(~behaviorP));
    
    
    
    
        rTau(i)=1;%min(20,find(diff(autocorr(rtemp-smooth(rtemp,100)',150))>0,1,'first'));

    gTau=fmax;%rTau(i);
    
    

    RGcorr(i)=corr2JN(gtemp,rtemp);
    
    gtemp=gtemp-mean(gtemp);
    nTerms=floor(length(gtemp)/gTau);
    [~, perm] = sort(rand(bootN,nTerms),2);

  %  perm=perm(:,rePerm);
    gperm=gtemp(perm);
    rperm=rtemp(perm);
    
   % gperm=randsample(gDistr,nTerms*bootN,1);
   % gperm=reshape(gperm,bootN,nTerms);
     rePerm=linspace(1,nTerms+.999,length(gtemp));
    rePerm=floor(rePerm);   
    gperm=gperm(:,rePerm);
 %   gperm=bsxfun(@minus,gperm, mean(gperm,2));
    gpermSTD=sqrt(mean(gperm.^2,2));
% behaviorF=behaviorF-mean(behaviorF);
% behaviorB=behaviorB-mean(behaviorB);
% behaviorP=behaviorP-mean(behaviorP);
% behaviorT=behaviorT-mean(behaviorT);
%rpermSTD=std(rtemp);
    GbootCorrF(i,:)=gperm*(behaviorF-mean(behaviorF))'./(std(behaviorF).*gpermSTD)/length(gtemp);
  %  RbootCorrF=corr2JN(rperm,behaviorF);
    GbootCorrB(i,:)=gperm*(behaviorB-mean(behaviorB))'./(std(behaviorB).*gpermSTD)/length(gtemp);
  %  RbootCorrB=corr2JN(rperm,behaviorB);
    GbootCorrP(i,:)=gperm*(behaviorP-mean(behaviorP))'./(std(behaviorP).*gpermSTD)/length(gtemp);
  %  RbootCorrP=corr2JN(rperm,behaviorP);
    GbootCorrT(i,:)=gperm*(behaviorT-mean(behaviorT))'./(std(behaviorT).*gpermSTD)/length(gtemp);
  %  RbootCorrT=corr2JN(rperm,behaviorT);
  
        deltabootGB(i,:)=mean((gperm(:,behaviorB)),2)-mean((gperm(:,~behaviorB)),2);
    deltabootGF(i,:)=mean((gperm(:,behaviorF)),2)-mean((gperm(:,~behaviorF)),2);
    deltabootGT(i,:)=mean((gperm(:,behaviorT)),2)-mean((gperm(:,~behaviorT)),2);
    deltabootGP(i,:)=mean((gperm(:,behaviorP)),2)-mean((gperm(:,~behaviorP)),2);

    

   % for iTrials=1:bootN
        
      %  gperm=gtemp(randperm(length(gtemp)));
%         plot(hasPointsTime,gperm+2);
%         hold on
%         area(hiResFrameTime,(ethogram==-1),'facecolor',bcolor,'edgecolor','none')
% hold on
% area(hiResFrameTime,(ethogram==0),'facecolor',pausecolor,'edgecolor','none')
% area(hiResFrameTime,(ethogram==1),'facecolor',fcolor,'edgecolor','none')
% area(hiResFrameTime,(ethogram==4),'facecolor',turncolor,'edgecolor','none')
% hold off
%         drawnow
    %    rperm=rtemp(randperm((length(gtemp))));
  % RGbootCorr(iTrials)=corr2JN(rperm,gperm);
    
%     deltabootRB(i,iTrials)=nanmean(rperm(behaviorB))-nanmean(rperm(~behaviorB));
%     deltabootGB(i,iTrials)=nanmean(gperm(behaviorB))-nanmean(gperm(~behaviorB));
%     deltabootRF(i,iTrials)=nanmean(rperm(behaviorF))-nanmean(rperm(~behaviorF));
%     deltabootGF(i,iTrials)=nanmean(gperm(behaviorF))-nanmean(gperm(~behaviorF));
%     deltabootRT(i,iTrials)=nanmean(rperm(behaviorT))-nanmean(rperm(~behaviorT));
%     deltabootGT(i,iTrials)=nanmean(gperm(behaviorT))-nanmean(gperm(~behaviorT));
%     
%        deltabootRP(i,iTrials)=nanmean(rperm(behaviorP))-nanmean(rperm(~behaviorP));
%     deltabootGP(i,iTrials)=nanmean(gperm(behaviorP))-nanmean(gperm(~behaviorP));
%     

    %end
%     GbootCorrSTDB(i)=quantile(GbootCorrB,.95);
%     RbootCorrSTDB(i)=quantile(RbootCorrB,.95);
%     GbootCorrSTDF(i)=quantile(GbootCorrF,.95);
%     RbootCorrSTDF(i)=quantile(RbootCorrF,.95);
%     GbootCorrSTDP(i)=quantile(GbootCorrP,.95);
%     RbootCorrSTDP(i)=quantile(RbootCorrP,.95);
%     GbootCorrSTDT(i)=quantile(GbootCorrT,.95);
%     RbootCorrSTDT(i)=quantile(RbootCorrT,.95);
    
%     GbootDeltaSTDB(i)=quantile(deltabootGB,.99);
%     RbootDeltaSTDB(i)=quantile(deltabootRB,.99);
%     GbootDeltaSTDF(i)=quantile(deltabootGF,.99);
%     RbootDeltaSTDF(i)=quantile(deltabootRF,.99);
%     GbootDeltaSTDT(i)=quantile(deltabootGT,.99);
%     RbootDeltaSTDT(i)=quantile(deltabootRT,.99);
%     
%     
  %  RGbootCorrSTD(i)=std(RGbootCorr);
    
    
end

%%
tic
bootN2=100000;
    gperm=randsample(gDistr,nTerms*bootN2,1);
    gperm=reshape(gperm,bootN2,nTerms);
     rePerm=linspace(1,nTerms+.999,length(gtemp));
    rePerm=floor(rePerm);   
    gperm=gperm(:,rePerm);
    gperm=bsxfun(@minus,gperm, mean(gperm,2));
    gpermSTD=sqrt(mean(gperm.^2,2));

    
    deltabootGB=mean((gperm(:,behaviorB)),2)-mean((gperm(:,~behaviorB)),2);
    deltabootGF=mean((gperm(:,behaviorF)),2)-mean((gperm(:,~behaviorF)),2);
    deltabootGT=mean((gperm(:,behaviorT)),2)-mean((gperm(:,~behaviorT)),2);
    deltabootGP=mean((gperm(:,behaviorP)),2)-mean((gperm(:,~behaviorP)),2);

    
    toc



%%

    GbootCorrSTDB=quantile(GbootCorrB,.995,2)';
    GbootCorrSTDF=quantile(GbootCorrF,.995,2)';
    GbootCorrSTDP=quantile(GbootCorrP,.995,2)';
    GbootCorrSTDT=quantile(GbootCorrT,.995,2)';

    
    deltaGBThresh=quantile(deltabootGB(:),deltaCI);
    deltaGFThresh=quantile(deltabootGF(:),deltaCI);
    deltaGTThresh=quantile(deltabootGT(:),deltaCI);
    deltaGPThresh=quantile(deltabootGP(:),deltaCI);
    s=7.5;
       
    deltaGBThresh=std(deltabootGB(:))*s;
    deltaGFThresh=std(deltabootGF(:))*s;
    deltaGTThresh=std(deltabootGT(:))*s;
    deltaGPThresh=std(deltabootGP(:))*s;
    [deltaGFThresh deltaGTThresh deltaGBThresh deltaGPThresh]
     
    
RGThresh=1;
deltaCI=.99999;
GRange=range(G2,2);
gthresh=.3;
maxN=17;
nanThresh=.1;
nancount=mean(isnan(G2),2)';
%possibleZ=find(abs(gzcorr)>.5);

possibleB=find(((GcorrB)>GbootCorrSTDB) & GcorrB>gthresh & RGcorr<RGThresh...
    & deltaGB>deltaGBThresh & nancount<nanThresh);
 [~,ia]=sort(GRange(possibleB),'descend');
 possibleB=possibleB(ia);
 possibleB=possibleB(1:min(length(possibleB),maxN));
possibleF=find(((GcorrF)>GbootCorrSTDF) & (GcorrF)>gthresh & RGcorr<RGThresh...
    & deltaGF>deltaGFThresh & nancount<nanThresh);
%possibleF=1:length(GcorrF);
[~,ia]=sort(GcorrF(possibleF),'descend');
 possibleF=possibleF(ia);
  possibleF=possibleF(1:min(length(possibleF),maxN));

possibleP=find(((GcorrP)>GbootCorrSTDP) & GcorrP>gthresh & RGcorr<RGThresh ...
    & deltaGP>deltaGPThresh & nancount<nanThresh);
 [~,ia]=sort(GRange(possibleP),'ascend');
 possibleP=possibleP(ia);
 possibleP=possibleP(max(length(possibleP)-maxN,1):end);
possibleB=possibleP;
 possibleT=find(((GcorrT)>GbootCorrSTDT) & GcorrT>gthresh & RGcorr<RGThresh...
    & deltaGT>deltaGTThresh & nancount<nanThresh);
 [~,ia]=sort(GRange(possibleT),'ascend');
 possibleT=possibleT(ia);
 possibleT=possibleT(max(length(possibleT)-maxN,1):end);
 
 overlap=intersect(possibleB,possibleT);
 possibleB=possibleB(~ismember(possibleB,overlap));
 possibleT=possibleT(~ismember(possibleT,overlap));
 moreB=deltaGP(overlap)>deltaGT(overlap);
possibleB=[possibleB overlap(moreB)];
possibleT=[possibleT overlap(~moreB)];

 
 %possibleB=possibleB(~ismember(possibleB,[possibleT,possibleF]));

 %possibleT=possibleT(~ismember(possibleT,[possibleB,possibleF]));
 
 %%
 
figure
spaceVec=[ 1 1 1];
i=-2;
plotOffset=0;
alphaVal=.2;
timeLength=max(hasPointsTime)+10;
ethoSelect=([1;diff(ethogram)]~=0 );
ethoSelect=ethoSelect(~isnan(ethogram));
FplotSelect=~isnan(ethoTrack);

ethoSelect=smooth(ethoSelect,3)>0;
ethoSelect(end)=true;
    ethoHeight=10;

for class=1:3
    space=spaceVec(class);
    plotOffset=i+2+plotOffset;
    plotOffset=-.5;
    switch class
        case 1
select=[possibleF ];
        case 2
            select=possibleT;
        case 3
            select=possibleB;
    end
    
    
%select=select(1:min(6,length(select)));
if ~isempty(select)
for i=1:length(select)
        subplot(3,2,2*(class-1)+1);
        plot(frameTimeTrack(FplotSelect),medfilt1(G2(select(i),FplotSelect),5)'+space*(i+plotOffset),'black','linew',1)

hold on

    subplot(3,2,2*(class-1)+2);
plot(frameTimeTrack(FplotSelect),R2(select(i),FplotSelect)'+space*(i+plotOffset),'black','linew',1)
hold on
%plot(R2(possibleB(i),:)'+space*i,'r')

end

    subplot(3,2,2*(class-1)+2);
ylim([0 ethoHeight])
set(gca,'YTick',space,'fontsize',12)
xlabel('Time(s)')
ylabel('\Delta F /F0')
  text( max(frameTimeTrack(FplotSelect)+5)*ones(size(select)), space*plotOffset+[space:space:space*(i)],cellstr(num2str(cgIdxRev(select))),'VerticalAlignment'...
        ,'middle', 'HorizontalAlignment','left','color',[0 0 0],...
        'fontsize',13);
    xlim([0 timeLength])

h(1)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)<0),'facecolor',bcolor);
h(2)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==0),'facecolor',pausecolor);
h(3)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==1),'facecolor',fcolor);
h(4)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==4),'facecolor',turncolor);
drawnow
for ii=1:length(h)
    
    subplot(3,2,2*(class-1)+1);
    h(ii).EdgeColor='none';
    h(ii).Face.ColorType='truecoloralpha';
    h(ii).Face.ColorData(4)=alphaVal*255;
end

ylim([0 ethoHeight])
set(gca,'YTick',space,'fontsize',12)
xlabel('Time(s)')
ylabel('\Delta F /F0')
  text( max(frameTimeTrack(FplotSelect)+5)*ones(size(select)), space*plotOffset+[space:space:space*(i)],cellstr(num2str(cgIdxRev(select))),'VerticalAlignment'...
        ,'middle', 'HorizontalAlignment','left','color',[0 0 0],...
        'fontsize',13);
    xlim([0 timeLength])
h(1)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)<0),'facecolor',bcolor);
h(2)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==0),'facecolor',pausecolor);
h(3)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==1),'facecolor',fcolor);
h(4)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==4),'facecolor',turncolor);
drawnow
for i=1:length(h)
    

    h(i).EdgeColor='none';
    h(i).Face.ColorType='truecoloralpha';
    h(i).Face.ColorData(4)=alphaVal*255;
end
end

end