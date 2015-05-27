control={'O:\20150311\BrainScanner20150311_143121\'};

gcamp={'F:\20141212\BrainScanner20141212_145951\',...
%'O:\20150118\BrainScanner20150118_184857 - 4+min\',...
'F:\20141214\BrainScanner20141214_175105 - Fred worm 2\'}
%%
corrPlotAll=[];RGall=[];
accumfunction = @(x) -qmean(-x,.2);
accumfunction = @(x) mean(x);
%mutualInformation=@(x) nansum(nansum(x.*log(x./(bsxfun(@times,nansum(x,1),nansum(x,2))))));

figure
for i=1:length(gcamp)
Gdata=load([gcamp{i} 'heatData'],'G2','R2','acorr','DmatAll','cgIdx','cgIdxRev');
acorr=Gdata.acorr;
DmatAll=Gdata.DmatAll;
rejects=isnan(Gdata.G2(:)+Gdata.R2(:));
G2temp=Gdata.G2(~rejects);
R2temp=Gdata.R2(~rejects);

%acorr=pdist2(Gdata.G2,Gdata.R2,@mutualInformation);
[corrRG]=histcn([G2temp,R2temp],-1:.1:2,-1:.1:2);

corrRG=corrRG/sum(corrRG(:));

MI=nansum(nansum(corrRG.*log2(corrRG./(bsxfun(@times,nansum(corrRG,1),nansum(corrRG,2))))));


G2temp=G2temp-mean(G2temp);
R2temp=R2temp-mean(R2temp);
corr2(G2temp,R2temp);
[coeff,score,latent,tsquared,explained,mu] = pca([G2temp R2temp]);

[MI,explained(2)]

%%
corrEdge=-1:.1:1;
dEdge=0:10:400;
DmatAll=DmatAll(:);
bad=DmatAll==0;
acorr=acorr(~bad);
DmatAll=DmatAll(~bad);

[corrD edges mid loc]=histcn([DmatAll(:),acorr(:)],dEdge,-1:.1:2);

corrPlot=accumarray(loc(:,1),acorr(:),size(mid{1}'),accumfunction);
corrPlotSTD=accumarray(loc(:,1),acorr(:),size(mid{1}'),@std);
corrPlotN=accumarray(loc(:,1),ones(size(acorr(:))),size(mid{1}'))/2;
corrPlotAll=[corrPlotAll,corrPlot];
plot(mid{1},corrPlot)
%errorbar(mid{1},corrPlot,corrPlotSTD./sqrt(corrPlotN),'b');
hold on
%%
try
cAll= load([gcamp{i}  filesep 'corrandPossibleCoor']);
c(i)=mean(cAll.RGcorr);
ci(i)=std(cAll.RGcorr);
cn(i)=length(cAll.RGcorr);
RGall=[RGall, cAll.RGcorr];
catch
    cAll=load([gcamp{i}  filesep 'RGcorr']);
    c(i)=mean(cAll.RGcorr);
ci(i)=std(cAll.RGcorr);
cn(i)=length(cAll.RGcorr);
RGall=[RGall, cAll.RGcorr];

end

end

%plot(mid{1},mean(corrPlotAll,2));
hold on

%%
for i2=1:length(control)
Gdata=load([control{i2} 'heatData'],'G2','R2','acorr','DmatAll','cgIdx','cgIdxRev');
acorr=Gdata.acorr;
DmatAll=Gdata.DmatAll;
rejects=isnan(Gdata.G2(:)+Gdata.R2(:));
G2temp=Gdata.G2(~rejects);
R2temp=Gdata.R2(~rejects);
acorr=pdist2(Gdata.G2,Gdata.R2,@mutualInformation);

[corrRG]=histcn([G2temp,R2temp],-1:.1:2,-1:.1:2);

corrRG=corrRG/sum(corrRG(:));
MI=nansum(nansum(corrRG.*log2(corrRG./(bsxfun(@times,nansum(corrRG,1),nansum(corrRG,2))))));

%mutualInformation(corrRG)
G2temp=G2temp-mean(G2temp);
R2temp=R2temp-mean(R2temp);
corr2(G2temp,R2temp);
[coeff,score,latent,tsquared,explained,mu] = pca([G2temp R2temp]);
[MI,explained(2)]
%%
corrEdge=-1:.1:1;
dEdge=0:10:400;
DmatAll=DmatAll(:);
bad=DmatAll==0;
acorr=acorr(~bad);
DmatAll=DmatAll(~bad);
[corrD edges mid loc]=histcn([DmatAll(:),acorr(:)],dEdge,-1:.1:1);

corrPlot=accumarray(loc(:,1),acorr(:),size(mid{1}'),accumfunction);
corrPlotSTD=accumarray(loc(:,1),acorr(:),size(mid{1}'),@std);
corrPlotN=accumarray(loc(:,1),ones(size(acorr(:))),size(mid{1}'))/2;

%errorbar(mid{1},corrPlot,corrPlotSTD./sqrt(corrPlotN));
plot(mid{1},corrPlot);

hold on
%%
try
cAll= load([control{i2}  filesep 'corrandPossibleCoor']);
c(i+i2)=mean(cAll.RGcorr);
ci(i+i2)=std(cAll.RGcorr);
cn(i+i2)=length(cAll.RGcorr);
catch
    cAll=load([control{i2}  filesep 'RGcorr']);
    c(i+i2)=mean(cAll.RGcorr);
ci(i+i2)=std(cAll.RGcorr);
cn(i+i2)=length(cAll.RGcorr);

end

end

%%
bar(1:4, c);
hold on
errorbar(1:4,c,ci./(sqrt(cn-1)),'.')

figure
plot(mid{1},mean(corrPlotAll,2));
hold on
plot(mid{1},corrPlot)

%%
figure
histc(RGall,-1:.2:1)
figure
histc(cAll.RGcorr,-1:.2:1)
