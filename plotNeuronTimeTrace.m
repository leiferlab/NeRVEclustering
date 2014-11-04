function plotNeuronTimeTrace(trackData,plotIdx,fluor) 
%function takes trackData from the runTrack program and plots time courses
%of the fluor for certain indices;
%plotIdx=[1,3,4,6,8,11,14,15,16,19,21,22,24,25,26,28,30,32,34];
figure
normalizeFlag=1;
interpFlag=1;
smoothWindow=3;
timeStep=.25;
if nargin==2
    fluor='red';
end
    switch fluor
        case 'green'
            output=trackData(:,5);
        case 'red'
            output=trackData(:,4);
        case 'ratio'
            output=trackData(:,4)./trackData(:,5);
        case 'x'
            output=trackData(:,1);
    end

    
    %for now, until i impliment output2
    output2=output;
    for i=1:length(plotIdx);
        idx=plotIdx(i);
        t=trackData((trackData(:,end)==idx),end-1);
        a=output((trackData(:,end)==idx));
        a2=output2((trackData(:,end)==idx));
        
   tnew=min(t):max(t);

        if interpFlag
          a=interp1(t,a,tnew);
          t=tnew;
          a=(smooth(a,smoothWindow));

        else
                   anew=nan(size(tnew));
            anew2=anew;
            anew(t-min(t)+1)=a;
     
         a=anew;
         anew2(t-min(t)+1)=a2;
         a2=anew2;

         t=tnew;
        end
        
        
        a0=nanmean(a(60:end));
        if normalizeFlag
     %        a=(smooth(a/nanmedian(a2),smoothWindow))/5+i-1;
    %         a2=(smooth(a2/nanmedian(a2),smoothWindow))/5+i-1;
a=((a-a0)/a0)+2*i;
%a=(a-median(a))/median(a);
        else
            a=(smooth(a,smoothWindow))+i-1;
            a2=(smooth(a2,smoothWindow))+i-1;
        end
        t=t*timeStep;
        plot(t,a,'black','linewidth',2);
        
        hold('on');
 %       plot(handles.axes2,t,a2,'g');
    end
xlabel('Time (s)');
xlim([0,max(t)]);
set(gca,'yticklabel',[]);
