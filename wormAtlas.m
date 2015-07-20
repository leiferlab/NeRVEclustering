atlas=load('Y:\Jeff\NeuronPositions');
backNeurons={'AVA','AVE','AIB','RIM','VA1'};
forwardNeurons={'AVB','VB1'};
deepVentralNeurons={'DD1','RIV','VA4','VB5','VC6','DA4','RIB','SMD'};


A=cellfun(@(x)strncmp(x,backNeurons,3),atlas.ID,...
    'uniformOutput',0);
A=cell2mat(A);
backIdx=any(A,2);

A=cellfun(@(x)  strncmp(x,forwardNeurons,3),atlas.ID,...
    'uniformOutput',0);
A=cell2mat(A);
forwardIdx=any(A,2);

A=cellfun(@(x)strncmp(x,deepVentralNeurons,3),atlas.ID,...
    'uniformOutput',0);
A=cell2mat(A);
deepIdx=any(A,2);


%%
figure
select=atlas.y<-2.3 & atlas.x<0.1;

colorBalls=.5+zeros(length(select),3);
colorBalls(deepIdx(select),3)=1;
colorBalls(backIdx(select),1)=1;
colorBalls(forwardIdx(select),2)=1;
transp=deepIdx(select)|backIdx(select)|forwardIdx(select);
transp=double(transp);
transp(transp==1)=.8;
transp(transp==0)=.2;
%subplot(2,1,1)
subplot(2,1,1);
scatter3sph(atlas.y(select)/7,atlas.z(select)/7,atlas.x(select)/7,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight

ax1=gca;
subplot(2,1,2)
select2=find(select);
select2=select2(transp==.8);

scatter3sph(atlas.y(select2)/7,atlas.z(select2)/7,atlas.x(select2)/7,...
    'size',.002,'color',colorBalls(transp==.8,:),'transp',.8);axis equal off tight
ax2=gca;
set([ax1 ax2],'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
    ,'manual','YLimMode','manual','ALimMode','manual');
 linkprop([ax1,ax2],{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle',...
     'Xlim','YLim','ZLim'})


% text(atlas.y(select2)/10,atlas.z(select2)/10,atlas.x(select2)/10,atlas.ID(select2),...
%     'HorizontalAlignment','left','FontSize',12);
% 

% subplot(2,1,2)
% 
% scatter3sph(fiducials(:,1),fiducials(:,2),3*fiducials(:,3),...
%     'size',10,'color',colorBalls,'transp',.65)
%  
% axis equal off
%% Fred's worm 20141214
dataFolder='F:\20141214\BrainScanner20141214_175105 - Fred worm 2';
masterData=load([dataFolder filesep 'refStackStraight']);
fiducials=masterData.masterFiducials;

data=load([dataFolder filesep 'heatData']);
possibleCorrData=load([dataFolder filesep 'corrandPossibleCoor2']);
possibleB=possibleCorrData.possibleB;
possibleT=possibleCorrData.possibleT;
possibleP=possibleCorrData.possibleP;
possibleF=possibleCorrData.possibleF;

%%
figure
% possibleB=[6    65    48    27    40    29];
% possibleF=[     45    46 ];
% possibleP=[  24    37    13    11    75    59];

colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleT,3)=1;
colorBalls(possibleB,1)=1;
colorBalls(possibleF,2)=1;

figure;subplot(2,1,1);
transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
% set(gca,'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
%     ,'manual','YLimMode','manual','ALimMode','manual');
ax1=gca;
subplot(2,1,2);
transp=.2*zeros(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
ax2=gca;
set([ax1 ax2],'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
    ,'manual','YLimMode','manual','ALimMode','manual');
 linkprop([ax1,ax2],{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle',...
     'Xlim','YLim','ZLim'})

%%
figure

[~,cgIdxRev]=sort(data.cgIdx);
cgIdxRev=cgIdxRev';
if ~isempty(possibleF)
scatter3sph(fiducials(possibleF,1)/4000,fiducials(possibleF,2)/4000,3*fiducials(possibleF,3)/4000,...
    'size',.002,'color',[0 1 0],'transp',.8);axis equal off tight
text(fiducials(possibleF,1)/4000+.004,fiducials(possibleF,2)/4000+.004,3*fiducials(possibleF,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleF),1,possibleF./possibleF)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
end
hold on
scatter3sph(fiducials(possibleB,1)/4000,fiducials(possibleB,2)/4000,3*fiducials(possibleB,3)/4000,...
    'size',.002,'color',[1 0 0],'transp',.8);axis equal off tight
scatter3sph(fiducials(possibleT,1)/4000,fiducials(possibleT,2)/4000,3*fiducials(possibleT,3)/4000,...
    'size',.002,'color',[0 0 1],'transp',.8);axis equal off tight
text(fiducials(possibleT,1)/4000+.004,fiducials(possibleT,2)/4000+.004,3*fiducials(possibleT,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleT),1,possibleT./possibleT)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
text(fiducials(possibleB,1)/4000+.004,fiducials(possibleB,2)/4000+.004,3*fiducials(possibleB,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleB),1,possibleB./possibleB)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);

text(atlas.y(select2)/10,atlas.z(select2)/10,atlas.x(select2)/10,atlas.ID(select2),...
    'HorizontalAlignment','left','FontSize',12);


%% make a movie
G2=data.G2;
G2ColorsLookup=G2;
G2ColorsLookup=G2ColorsLookup/1.25;
G2ColorsLookup(G2ColorsLookup>1)=1;
G2ColorsLookup=colNanFill(G2ColorsLookup')';
G2ColorsLookup=ceil(G2ColorsLookup*100);

G2ColorsLookup(isnan(G2ColorsLookup))=1;
G2ColorsLookup(G2ColorsLookup<1)=1;

Pmap=cbrewer('seq','Blues',100);
Bmap=cbrewer('seq','YlOrRd',100);
Fmap=cbrewer('seq','Greens',100);

for t=1:size(G2,2);
   colorT=G2ColorsLookup(:,t);

   colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleP,:)=Pmap(G2ColorsLookup(possibleP,t),:);
colorBalls(possibleB,:)=Pmap(G2ColorsLookup(possibleB,t),:);
colorBalls(possibleF,:)=Pmap(G2ColorsLookup(possibleF,t),:);

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleP])=.8;

delete(gca)
scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
pause(.02)
end

%% Fred's worm 20150118 LONG dataset

dataFolder='O:\20150118\BrainScanner20150118_184857 - 4+min' ;
masterData=load([dataFolder filesep 'refStackStraight']);

masterData=load([dataFolder filesep 'refStackStraight']);
fiducials=masterData.masterFiducials;

data=load([dataFolder filesep 'heatData']);
possibleCorrData=load([dataFolder filesep 'corrandPossibleCoor2']);
possibleB=possibleCorrData.possibleB;
possibleT=possibleCorrData.possibleT;
possibleP=possibleCorrData.possibleP;
possibleF=possibleCorrData.possibleF;

%%
figure
% possibleB=[6    65    48    27    40    29];
% possibleF=[     45    46 ];
% possibleP=[  24    37    13    11    75    59];

colorBalls=.3+zeros(length(fiducials),3);
colorBalls(possibleT,3)=1;
colorBalls(possibleB,1)=1;
colorBalls(possibleF,2)=1;

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

%%
figure;subplot(2,1,1);
transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,30*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
% set(gca,'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
%     ,'manual','YLimMode','manual','ALimMode','manual');
ax1=gca;
subplot(2,1,2);
transp=.2*zeros(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,30*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
ax2=gca;
set([ax1 ax2],'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
    ,'manual','YLimMode','manual','ALimMode','manual');
 linkprop([ax1,ax2],{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle',...
     'Xlim','YLim','ZLim'})

%%
figure

[~,cgIdxRev]=sort(data.cgIdx);
cgIdxRev=cgIdxRev';
scatter3sph(fiducials(possibleF,1)/4000,fiducials(possibleF,2)/4000,30*fiducials(possibleF,3)/4000,...
    'size',.002,'color',[0 1 0],'transp',.8);axis equal off tight
hold on
scatter3sph(fiducials(possibleB,1)/4000,fiducials(possibleB,2)/4000,30*fiducials(possibleB,3)/4000,...
    'size',.002,'color',[1 0 0],'transp',.8);axis equal off tight
scatter3sph(fiducials(possibleT,1)/4000,fiducials(possibleT,2)/4000,30*fiducials(possibleT,3)/4000,...
    'size',.002,'color',[0 0 1],'transp',.8);axis equal off tight
text(fiducials(possibleF,1)/4000+.004,fiducials(possibleF,2)/4000+.004,3*fiducials(possibleF,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleF),1,possibleF./possibleF)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
text(fiducials(possibleT,1)/4000+.004,fiducials(possibleT,2)/4000+.004,3*fiducials(possibleT,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleT),1,possibleT./possibleT)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
text(fiducials(possibleB,1)/4000+.004,fiducials(possibleB,2)/4000+.004,3*fiducials(possibleB,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleB),1,possibleB./possibleB)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);

% text(atlas.y(select2)/10,atlas.z(select2)/10,atlas.x(select2)/10,atlas.ID(select2),...
%     'HorizontalAlignment','left','FontSize',12);


%% make a movie
G2=data.G2;
G2ColorsLookup=G2;
G2ColorsLookup=G2ColorsLookup/1.25;
G2ColorsLookup(G2ColorsLookup>1)=1;
G2ColorsLookup=colNanFill(G2ColorsLookup')';
G2ColorsLookup=ceil(G2ColorsLookup*100);

G2ColorsLookup(isnan(G2ColorsLookup))=1;
G2ColorsLookup(G2ColorsLookup<1)=1;

Pmap=cbrewer('seq','Blues',100);
Bmap=cbrewer('seq','YlOrRd',100);
Fmap=cbrewer('seq','Greens',100);

for t=1:size(G2,2);
   colorT=G2ColorsLookup(:,t);

   colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleP,:)=Pmap(G2ColorsLookup(possibleP,t),:);
colorBalls(possibleB,:)=Pmap(G2ColorsLookup(possibleB,t),:);
colorBalls(possibleF,:)=Pmap(G2ColorsLookup(possibleF,t),:);

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleP])=.8;

delete(gca)
scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
pause(.02)
end


%% Jeff's worm 20141212

dataFolder='F:\20141212\BrainScanner20141212_145951' ;

masterData=load([dataFolder filesep 'refStackStraight']);
fiducials=masterData.masterFiducials;

data=load([dataFolder filesep 'heatData']);
possibleCorrData=load([dataFolder filesep 'corrandPossibleCoor2']);
possibleB=possibleCorrData.possibleP;
possibleT=possibleCorrData.possibleT;
possibleP=possibleCorrData.possibleP;
possibleF=possibleCorrData.possibleF;


%%
figure
% possibleB=[2    61    67];
% possibleF=    [20    18    60    10] ;%old result:[20    18    42    10];
% possibleP=[14    19    35    29    33    45     7];%[  19    35    29    64    33    45     7];

colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleT,3)=1;
colorBalls(possibleB,1)=1;

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;


colorBalls(possibleF,2)=1;
%%
figure;subplot(2,1,1);
transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
% set(gca,'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
%     ,'manual','YLimMode','manual','ALimMode','manual');
ax1=gca;
subplot(2,1,2);
transp=.2*zeros(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
ax2=gca;
set([ax1 ax2],'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
    ,'manual','YLimMode','manual','ALimMode','manual');
 linkprop([ax1,ax2],{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle',...
     'Xlim','YLim','ZLim'})

%%
figure
[~,cgIdxRev]=sort(data.cgIdx);
cgIdxRev=cgIdxRev';
scatter3sph(fiducials(possibleF,1)/4000,fiducials(possibleF,2)/4000,3*fiducials(possibleF,3)/4000,...
    'size',.002,'color',[.3 1 .3],'transp',.8);axis equal off tight
hold on

scatter3sph(fiducials(possibleB,1)/4000,fiducials(possibleB,2)/4000,3*fiducials(possibleB,3)/4000,...
    'size',.002,'color',[1 .3 .3],'transp',.8);axis equal off tight
scatter3sph(fiducials(possibleT,1)/4000,fiducials(possibleT,2)/4000,3*fiducials(possibleT,3)/4000,...
    'size',.002,'color',[.3 .3 1],'transp',.8);axis equal off tight
text(fiducials(possibleF,1)/4000+.004,fiducials(possibleF,2)/4000+.004,3*fiducials(possibleF,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleF),1,possibleF./possibleF)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
text(fiducials(possibleT,1)/4000+.004,fiducials(possibleT,2)/4000+.004,3*fiducials(possibleT,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleT),1,possibleT./possibleT)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
text(fiducials(possibleB,1)/4000+.004,fiducials(possibleB,2)/4000+.004,3*fiducials(possibleB,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleB),1,possibleB./possibleB)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);



%% GFP control worm. 

dataFolder='O:\20150311\BrainScanner20150311_143121' ;

masterData=load([dataFolder filesep 'refStackStraight']);
fiducials=masterData.masterFiducials;
fiducials=fiducials(:,[1:2 4]);
data=load([dataFolder filesep 'heatData']);
possibleCorrData=load([dataFolder filesep 'corrandPossibleCoor2']);
possibleB=possibleCorrData.possibleB;
possibleT=possibleCorrData.possibleT;
possibleP=possibleCorrData.possibleP;
possibleF=possibleCorrData.possibleF;



figure
% possibleB=[2    61    67];
% possibleF=    [20    18    60    10] ;%old result:[20    18    42    10];
% possibleP=[14    19    35    29    33    45     7];%[  19    35    29    64    33    45     7];

colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleT,3)=1;
colorBalls(possibleB,1)=1;
colorBalls(possibleF,2)=1;

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;


figure;subplot(2,1,1);
transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,5*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
% set(gca,'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
%     ,'manual','YLimMode','manual','ALimMode','manual');
axis equal

ax1=gca;
subplot(2,1,2);
transp=.2*zeros(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,5*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
axis equal

ax2=gca;
set([ax1 ax2],'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
    ,'manual','YLimMode','manual','ALimMode','manual');
 linkprop([ax1,ax2],{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle',...
     'Xlim','YLim','ZLim'})


%% make a movie
G2=data.G2;
G2ColorsLookup=G2;
G2ColorsLookup=G2ColorsLookup/1.25;
G2ColorsLookup(G2ColorsLookup>1)=1;
G2ColorsLookup=colNanFill(G2ColorsLookup')';
G2ColorsLookup=bsxfun(@rdivide,G2ColorsLookup,max(G2ColorsLookup,[],2));

G2ColorsLookup=ceil(G2ColorsLookup*100);

G2ColorsLookup(isnan(G2ColorsLookup))=1;
G2ColorsLookup(G2ColorsLookup<1)=1;
Pmap=cbrewer('seq','YlGnBu',100);
Bmap=cbrewer('seq','YlOrRd',100);
Fmap=cbrewer('seq','YlGn',100);

for t=1:size(G2,2);
   colorT=G2ColorsLookup(:,t);

   colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleP,:)=Bmap(G2ColorsLookup(possibleP,t),:);
colorBalls(possibleB,:)=Bmap(G2ColorsLookup(possibleB,t),:);
colorBalls(possibleF,:)=Bmap(G2ColorsLookup(possibleF,t),:);

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleP])=.8;

delete(gca)
scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
view([90,90])

pause(.02)
end




%% anyFolder
dataFolder=uipickfiles('filterspec', 'V:');
dataFolder=dataFolder{1};
%%
masterData=load([dataFolder filesep 'tempFiducials']);
corrandPossibleC=uipickfiles('filterspec',dataFolder);

%%
fiducialsAll=masterData.fiducialPoints;

data=load([dataFolder filesep 'heatData']);
possibleCorrData=load(corrandPossibleC{1});
possibleB=possibleCorrData.possibleB;
possibleT=possibleCorrData.possibleT;
possibleP=possibleCorrData.possibleP;
possibleF=possibleCorrData.possibleF;

%%
masterIdx=80;
figur
fiducials=fiducialsAll{masterIdx};
fiducials(rejects(1:length(fiducials)),:)=[];
fiducials(cellfun(@(x) isempty(x),fiducials))={nan};
fiducials=cell2mat(fiducials);
fiducials=fiducials(:,[1 2 3]);
% possibleB=[6    65    48    27    40    29];
% possibleF=[     45    46 ];
% possibleP=[  24    37    13    11    75    59];

colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleT,3)=1;
colorBalls(possibleB,1)=1;
colorBalls(possibleF,2)=1;
fiducials([possibleB possibleT possibleF])

fiducials(:,2:3)=(-fiducials(:,2:3));
%%
figure;
transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

zFactor=50;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,zFactor*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
% set(gca,'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
%     ,'manual','YLimMode','manual','ALimMode','manual');
ax1=gca;
savefig([dataFolder filesep 'wormBallsAll'])
figure
transp=.2*zeros(1,length(fiducials));
transp([possibleB possibleF possibleT])=.8;

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,zFactor*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight
ax2=gca;
set([ax1 ax2],'cameraviewanglemode','manual','ZLimMode','manual','XLimMode'...
    ,'manual','YLimMode','manual','ALimMode','manual');
%  linkprop([ax1,ax2],{'CameraPosition','CameraUpVector','CameraTarget','CameraViewAngle',...
%      'Xlim','YLim','ZLim'})
%  
text(fiducials(possibleT,1)/4000+.004,fiducials(possibleT,2)/4000+.004,zFactor*fiducials(possibleT,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleT)',1,possibleT./possibleT)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);

text(fiducials(possibleB,1)/4000+.004,fiducials(possibleB,2)/4000+.004,zFactor*fiducials(possibleB,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleB)',1,possibleB./possibleB)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
text(fiducials(possibleF,1)/4000+.004,fiducials(possibleF,2)/4000+.004,zFactor*fiducials(possibleF,3)/4000+.004,...
    cellfun(@(x) num2str(x),(mat2cell(cgIdxRev(possibleF)',1,possibleF./possibleF)),'uniformOutput',0)...
    , 'HorizontalAlignment','left','FontSize',12);
savefig([dataFolder filesep 'wormBallsLabel'])



 
 