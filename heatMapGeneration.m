
function heatMapGeneration(dataFolder,RvalAll,GvalAll)


% take a datafolder, loading the wormFiducialIntensities to produce
% heatmaps and ordering for Red, green and ratios
if nargin==0
    dataFolder=uipickfiles();
    dataFolder=dataFolder{1};
end


%% fitting parameters 
%exponential fit equation
Fexponent=fittype('a*exp(b*x)+c','dependent',{'y'},'independent',...
    {'x'},'coefficients',{'a', 'b', 'c'});

fitOptions=fitoptions(Fexponent);
fitOptions.Lower=[0,-.2,0];
fitOptions.Upper=[1000,0,10000];

minWindow=150;
min_quant=30;


%% PHOTOBLEACHING CORRECTION

% intialize photobleaching corrections
photoBleachingR=zeros(size(RvalAll));
photoBleachingG=zeros(size(RvalAll));

% do photobleaching correction, correction is done by first running each R
% and G trace through an ordfilt function, filtering values to a 30th value
% out of 150 value window. 
for i=1:size(RvalAll,1)
    try
        %%
        %initialize x values for fitting y=f(x)
        xVals=(1:size(RvalAll,2))';
        % only take values where bot R and G are present
        present=(~isnan(RvalAll(i,:)+GvalAll(i,:))') ;
        xVals=xVals(present);
        % get R and G traces
        rVals=RvalAll(i,:)';
        gVals=GvalAll(i,:)';
        gVals=gVals(present);
        rVals=rVals(present);
        % do ord filtering 
        gVals=ordfilt2(gVals,min_quant,true(minWindow,1));
        rVals=ordfilt2(rVals,min_quant,true(minWindow,1));
        
        %set up more fitting parameters for Red, and fit starting point
        fitOptions.StartPoint=[range(rVals(rVals~=0)),-.0006,min(rVals(rVals~=0))];
        fitOptions.Weights=zeros(size(rVals));
        fitOptions.Weights(minWindow:end-minWindow)=1;
        
        %do exponential fitting
        [f,fout]=fit(xVals,rVals,Fexponent,fitOptions);
        
        %if fit is bad, try fit linear to loglinear plot
        if fout.rsquare<.9
            logVals=log(rVals);
            logVals=logVals(rVals~=0);
            logXvals=xVals(rVals~=0); %not actually logging xvals
            expFit=polyfit(logXvals,logVals,1);            
            f.a=exp(expFit(2));
            f.b=expFit(1);
        end
        
        %do the same for the green
        fitOptions.StartPoint=[range(gVals),-.001,min(gVals)];
        fitOptions.Weights=zeros(size(gVals));
        fitOptions.Weights(minWindow:end-minWindow)=1;

        %green always has a strange bump in intensity at the start, fit the
        %exponential starting after this by setting weights for the first
        %part to zero.
        [~,maxPos]=max(gVals(1:300));
        fitOptions.Weights(1:maxPos)=0;
        
        [g,gout]=fit(xVals,gVals,Fexponent,fitOptions);
        
        if f(1)>(max(RvalAll(i,:))+100)
            f=fit(xVals,rVals,'poly1');
            if f.p1>0
                f.p1=0;
            end
        end
        if g(1)>(max(GvalAll(i,:))+1000)
            g=fit(xVals,gVals,'poly1');
            if g.p1>0
                g.p1=0;
            end
        end
        %plot some of the results, turned off for now
        if 0 
            subplot(2,1,1);
            plot(GvalAll(i,:))
            hold on
            plot(g)
            ylim([0 g(0)+100])
            
            hold off
            subplot(2,1,2);
            plot(RvalAll(i,:))
            hold on
            
            plot(f)
            ylim([0 f(0)+100]);
            hold off
            drawnow
            pause(.1)
        end
        %calculating photobleaching correction from exponential fits
        photoBleachingR(i,:)=f((1:size(RvalAll,2)))-f(size(RvalAll,2));
        photoBleachingG(i,:)=g((1:size(RvalAll,2)))-g(size(RvalAll,2));
    catch me
        me
    end
    
    
end
%%
%apply photobleaching correction, nan the values that are very bright or
%dark 
rPhotoCorr=RvalAll-photoBleachingR ;
RvalstempZ=bsxfun(@minus,rPhotoCorr,nanmean(rPhotoCorr,2));
RvalstempZ=bsxfun(@rdivide,RvalstempZ,nanstd(RvalstempZ,[],2));
rPhotoCorr(RvalstempZ<-2|RvalstempZ>5|rPhotoCorr<40)=nan;


gPhotoCorr=GvalAll-photoBleachingG ;
GvalstempZ=bsxfun(@minus,gPhotoCorr,nanmean(gPhotoCorr,2));
GvalstempZ=bsxfun(@rdivide,GvalstempZ,nanstd(GvalstempZ,[],2));
gPhotoCorr(GvalstempZ>5|gPhotoCorr<0)=nan;


%% apply smoothing and fold change over baseline calculation

%Apply it to red
A=rPhotoCorr';
Asmooth=smooth2a(A,50,0);
Asmooth=colNanFill(Asmooth);
A0=quantile(Asmooth,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=imfilter(A2, gausswin(5,1)/sum( gausswin(5,1)));
A2(A2<-1)=-nan;
R2=A2';

%then to green
A=(gPhotoCorr)';
Asmooth=smooth2a(A,50,0);
Asmooth=colNanFill(Asmooth);
A0=quantile(Asmooth,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=imfilter(A2, gausswin(5,1)/sum( gausswin(5,1)));
A2(A2<-1)=-nan;
G2=A2';

%chop out flashes or other strange values

rmean=nanmean(R2(:));
rstd=nanstd(R2(:));
nanmapr=R2>4|isnan(R2);%(rmean+3*rstd);
gmean=nanmean(G2(:));
gstd=nanstd(G2(:));
nanmapg=G2>4|isnan(G2);%(gmean+3*gstd);
G2(nanmapg)=nan;
gPhotoCorr(nanmapg)=nan;
rPhotoCorr(nanmapr)=nan;

%now do it for the ratio of  R to G
gfilt=@(x,h) imfilter(x, gausswin(h,1)/sum( gausswin(h,1)));

%fill in nans, smooth both R and G, then take ratio
Gsmooth=colNanFill(gPhotoCorr');
Rsmooth=colNanFill(rPhotoCorr');
Gsmooth=gfilt(Gsmooth,5);
Rsmooth=gfilt(Rsmooth,5);
A=Gsmooth./Rsmooth;
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
Ratio2=A2';

%reinsert nans into Ratio, filling in some of the isolated nans. 
nan_map=isnan(gPhotoCorr+rPhotoCorr);
bad_col=mean(nan_map)>.3;
nan_map(:,bad_col)=1;
nan_map=imopen(nan_map,[ 1 1 1 ]);
nan_map(:,bad_col)=1;
nan_map=imclose(nan_map,ones(1,10));
Ratio2(nan_map)=nan;

%% sort rows of correlation matrix  using heirarchical clustering,

A(isnan(A))=0;
acorr=corr(A);
atemp=nancov(A)./sqrt(nanvar(A)'*nanvar(A));
acorr(isnan(acorr))=atemp(isnan(acorr));
acorr(isnan(acorr))=0;

cg = clustergram(acorr);
cgIdx=str2double(get(cg,'RowLabels'));
[~,cgIdxRev]=sort(cgIdx);
%close annoying clustergram plot
close all hidden

rRaw=RvalAll;
gRaw=GvalAll;


%%
save([dataFolder filesep 'heatData'],'G2','R2','gRaw','rRaw',...
    'rPhotoCorr','gPhotoCorr','Ratio2','acorr','cgIdx','cgIdxRev');

