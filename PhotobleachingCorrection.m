function [red_corr,green_corr]=PhotobleachingCorrection(red_signal,green_signal);

minWindow=150;
slidingWindow=30;
        
%intialize photobleaching traces
photoBleachingR=zeros(size(red_signal));
photoBleachingG=photoBleachingR;

%photobleaching model function
Fexponent=fittype('a*exp(b*x)+c','dependent',{'y'},'independent',...
    {'x'},'coefficients',{'a', 'b', 'c'});s

fitOptions=fitoptions(Fexponent);
fitOptions.Lower=[0,-.2,0];
fitOptions.Upper=[1000,0,10000];

%progressbar(0)
for i=1:size(red_signal,1)
    
    try
        %%
        
        % initialize x trace, remove nans
        xVals=(1:size(red_signal,2))';
        present=(~isnan(red_signal(i,:)+green_signal(i,:))') ;
        xVals=xVals(present);
        
        rVals=red_signal(i,:)';
        gVals=green_signal(i,:)';
        gVals=gVals(present);
        rVals=rVals(present);
        % sliding min filter
        gVals=ordfilt2(gVals,slidingWindow,true(minWindow,1));
        rVals=ordfilt2(rVals,slidingWindow,true(minWindow,1));
        
        %fit red trace
        fitOptions.StartPoint=[range(rVals(rVals~=0)),-.0006,min(rVals(rVals~=0))];
        fitOptions.Weights=zeros(size(rVals));
        fitOptions.Weights(minWindow:end-minWindow)=1;
        [f,fout]=fit(xVals,rVals,Fexponent,fitOptions);
        
        % if rsquare is bad, just fit a line to a log plot
        if fout.rsquare<.9
            logVals=log(rVals);
            logVals=logVals(rVals~=0);
            logXvals=xVals(rVals~=0);
            expFit=polyfit(logXvals,logVals,1);
            [f.b expFit(1) 1-fout.rsquare]
            
            f.a=exp(expFit(2));
            f.b=expFit(1);
            
        end
        
        %fit green trace
        fitOptions.StartPoint=[range(gVals),-.001,min(gVals)];
        fitOptions.Weights=zeros(size(gVals));
        fitOptions.Weights(minWindow:end-minWindow)=1;
        %unweight points before the initial peak, not sure why that peak is
        %there but it biases the fit.
        [~,maxPos]=max(gVals(1:300));
        fitOptions.Weights(1:maxPos)=0;
        [g,gout]=fit(xVals,gVals,Fexponent,fitOptions);

        % plot traces
        if 0
            subplot(2,1,1);
            plot(green_signal(i,:))
            hold on
            plot(g)
            ylim([0 g(0)+100])
            
          
            hold off
            subplot(2,1,2);
            plot(red_signal(i,:))
            hold on
            
            plot(f)
            ylim([0 f(0)+100]);
            hold off
            drawnow
            pause(.1)
        end
        photoBleachingR(i,:)=f((1:size(red_signal,2)))-f(size(red_signal,2));
        photoBleachingG(i,:)=g((1:size(red_signal,2)))-g(size(red_signal,2));
    catch me
        me
    end

end

red_corr=red_signal-photoBleachingR;
green_corr=green_signal-photoBleachingG;