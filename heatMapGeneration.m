
function heatMapGeneration(dataFolder,rRaw,gRaw)

%heatMapGeneration takes takes raw Red and Green signals from whole brain
%imaging and does processing to create ratiometric signals


%%% INPUTS
%dataFolder - destination folder for saving nerual signal results
%rRaw - an N by T matrix of red fluor signals for N neurons and T times
%gRaw - same as above but for the green signal.


%%%OUTPUTS
% inside dataFolder, a heatdata.mat file will be saved with signal at
% various parts of the processing. 

%	rRaw - an N neurons x T volumes matrix with the raw red signal from each
%of the neurons. Signals are averaged pixel values around each tracked 
%neuron with no other processing except for flash removal. 

%	gRaw - same as rRaw but for the green signal. 

%	rPhotoCorr - the rRaw signal after photobleaching correction for each 
%neuron. No other smoothing or normalization is applied. Photobleaching 
%correction is applied by fitting an exponential curve to a wide 20th 
%percentile filter, and then subtraction the exponential from the raw signal.

%	gPhotoCorr - same as above but with the green signal. 
%Exponential curves are fit independently. 

%	R2 - Smoothed and normalized version of rPhotoCorr. Normalization is %
%done as delta F/ F0, where F0 is the lower 20th percentile signal. 
%	G2 - Same as above but with gPhotoCorr.

%	Ratio2 - The ratio signal is defined as gPhotoCorr/rPhotoCorr, 
%the Ratio is then normalized as delta R/ R0. is the same way as R2 and G2. 

% cgIdx - indices used to organize output, obtained by doing heirarchical
% clustering of the correlation matrix of Ratio2.
% cgIdxRev - used to go back from organized output to original one. Not
% really necessary anymore. 

%In our current setup, it works out to a fwhm of about 5 steps, or slightly
%less than a second



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
photoBleachingR=zeros(size(rRaw));
photoBleachingG=zeros(size(gRaw));

% do photobleaching correction, correction is done by first running each R
% and G trace through an ordfilt function, filtering values to a 30th value
% out of 150 value window. 
for i=1:size(rRaw,1)
    try
        %%
        %initialize x values for fitting y=f(x)
        xVals=(1:size(rRaw,2))';
        % only take values where bot R and G are present
        present=(~isnan(rRaw(i,:)+gRaw(i,:))') ;
        present=present & (rRaw(i,:)~=0)' & (gRaw(i,:) ~=0)';
        xVals=xVals(present);
        % get R and G traces
        rVals=rRaw(i,:)';
        gVals=gRaw(i,:)';
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
        
        if f(1)>(max(rRaw(i,:))+100)
            f=fit(xVals,rVals,'poly1');
            if f.p1>0
                f.p1=0;
            end
        end
        if g(1)>(max(gRaw(i,:))+1000)
            g=fit(xVals,gVals,'poly1');
            if g.p1>0
                g.p1=0;
            end
        end
        %plot some of the results, turned off for now
        if 0
            subplot(2,1,1);
            plot(gRaw(i,:))
            hold on
            plot(g)
            ylim([0 g(0)+100])
            
            hold off
            subplot(2,1,2);
            plot(rRaw(i,:))
            hold on
            
            plot(f)
            ylim([0 f(0)+100]);
            hold off
            drawnow
            pause(.1)
        end
        
        limit=min(3000,size(rRaw,2));
        %calculating photobleaching correction from exponential fits
        photoBleachingR(i,:)=f((1:size(rRaw,2)))-f(limit);
        photoBleachingG(i,:)=g((1:size(rRaw,2)))-g(limit);
    catch me
        me
    end
    
    
end
%%
%apply photobleaching correction, nan the values that are very bright or
%dark 
rPhotoCorr=rRaw-photoBleachingR ;
RvalstempZ=bsxfun(@minus,rPhotoCorr,nanmean(rPhotoCorr,2));
RvalstempZ=bsxfun(@rdivide,RvalstempZ,nanstd(RvalstempZ,[],2));
rPhotoCorr(RvalstempZ<-2|RvalstempZ>5|rPhotoCorr<40)=nan;


gPhotoCorr=gRaw-photoBleachingG ;
GvalstempZ=bsxfun(@minus,gPhotoCorr,nanmean(gPhotoCorr,2));
GvalstempZ=bsxfun(@rdivide,GvalstempZ,nanstd(GvalstempZ,[],2));
gPhotoCorr(GvalstempZ>5|gPhotoCorr<0)=nan;


%% apply smoothing and fold change over baseline calculation

%Process red and green signals, functions below. 
R2=processSignal(rPhotoCorr);
G2=processSignal(gPhotoCorr);

%chop out flashes or other strange values
nanmapr=R2>4|isnan(R2);
nanmapg=G2>4|isnan(G2);
G2(nanmapg)=nan;
gPhotoCorr(nanmapg)=nan;
rPhotoCorr(nanmapr)=nan;

%now, process the ratio
Ratio2=processRatio(rPhotoCorr,gPhotoCorr);



%% sort rows of correlation matrix  using heirarchical clustering,
A=Ratio2';
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



%%
save([dataFolder filesep 'heatData'],'G2','R2','gRaw','rRaw',...
    'rPhotoCorr','gPhotoCorr','Ratio2','acorr','cgIdx','cgIdxRev');


function A_out=processSignal(A)

%do smoothing and deltaF/F0 processing for red and green signals

%define a filtering function. h will be the fwhm of the gaussian
gfilt=@(x,h) imfilter(x, gausswin(15,15/h)/sum( gausswin(15,15/h)));

%start with transpose for colNanFill to work
A=A';

%for finding the baseline, smooth and then take 20% quantile, this is not
%the smoothing used for the output.
Asmooth=smooth2a(A,50,0);
Asmooth=colNanFill(Asmooth);
A0=quantile(Asmooth,.2,1);

%take deltaF/F0
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);

%fill in nans and smooth
A2=colNanFill(A);
A2=gfilt(A2,5);
%very negative signals are likely bad
A2(A2<-1)=-nan;

A_out=A2';


function Aout=processRatio(R,G)
%same as above, but for Ratio. 
gfilt=@(x,h) imfilter(x, gausswin(15,15/h)/sum( gausswin(15,15/h)));

%fill in nans, smooth both R and G, then take ratio
Gsmooth=colNanFill(G');
Rsmooth=colNanFill(R');

Gsmooth=gfilt(Gsmooth,5);
Rsmooth=gfilt(Rsmooth,5);

A=Gsmooth./Rsmooth;

%find lower percentile for deltaR/R0
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);

A2=colNanFill(A);

Aout=A2';

%reinsert nans into Ratio, filling in some of the isolated nans. 

%find nans
nan_map=isnan(G+R);

%if more than .3 of the data in a col are nan, trash the col
bad_col=mean(nan_map)>.5;
nan_map(:,bad_col)=1;

%do morphological open, removing isolated nans, I'm ok 
%interpolating through some of these
nan_map=imopen(nan_map,[ 1 1 1 ]);
nan_map(:,bad_col)=1;

%if man nans appear, merge them. 
nan_map=imclose(nan_map,ones(1,10));

Aout(nan_map)=nan;
