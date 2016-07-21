%select file name

movieFile=uipickfiles('FilterSpec','Y:\PanNeuronal\20140819');

%% load vidObj using VideoReader
vidObj = VideoReader(movieFile{1});
lastFrame = read(vidObj, inf);
 numFrames= vidObj.NumberOfFrames;

%% load some number of stacks
clear imStack
clear unFilteredStack
stackCounter=1;
for iFrame=1:1:min(40,numFrames)
    
lastFrame = read(vidObj, iFrame);


lastFrame=normalizeRange(sum(double(lastFrame),3));
lastFramebp=bpass(lastFrame,.5,[30,30]);

imStack(:,:,stackCounter)=lastFramebp;
unFilteredStack(:,:,stackCounter)=lastFrame;
stackCounter=stackCounter+1;
end
%%
imagesc(max(imStack,[],3));
rect=round(getrect);
imStackCrop=imStack(rect(1):rect(1)+rect(3),rect(2):rect(2)+rect(4),:);

imStackSmooth=smooth3(imStackCrop,'box',[3,3,5]);

%%
smooth_weight = 1; 
image_weight = 100; 
delta_t = 4; 
%%
margin = 10; 
% phi = zeros(size(V)); 
% phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1; 
% phi = ac_reinit(phi-.5); 
% 
% 

% %%
% phi2=repmat(phi(:,:,20),[1,1,1])-max(phi(:))/2;
phi2=normalizeRange(V)-.5;
phi2=phi2(:,:,1);
%% testing single 2d plane



for iSlice=1:size(V,3)
    hold off;
    Vin=V(:,:,iSlice);
    Vin=smooth2a(Vin, 1,1);
   imagesc(Vin);
   
    hold on
    if iSlice==1;
        nIterations=20;
    else
        nIterations=10;
    end
    
for i = 1:1
    
    phi2 = ac_ChanVese_model(max(phi2(:))*(Vin), phi2, smooth_weight, image_weight, delta_t, nIterations); 
    
    if exist('h','var') && all(ishandle(h)),delete(h); end
  % h = patch(iso,'facecolor','w');  axis equal;  view(3); 
  [C,h]=contour(phi2,[0,0],'r');axis equal;
  
    set(gcf,'name', sprintf('#iters = %d',i));
    drawnow; 
end
contours(iSlice).C=C;
end



%% testing single3d plane
smooth_weight = .01; 
image_weight = 1000; 
delta_t = 4; 

margin = 10; 
% phi = zeros(size(V)); 
% phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1; 
% phi = ac_reinit(phi-.5); 
% 
gaussKernal=gausswin(30);
gaussKernal=convnfft(gaussKernal,gaussKernal');
gaussKernal=convnfft(gaussKernal,permute(gausswin(16),[2,3,1]));
gaussKernal=gaussKernal/sum(gaussKernal(:));
gaussKernal=gaussKernal-mean(gaussKernal(:));
worm2Smooth=imfilter(worm2,gaussKernal);


V2=normalizeRange(worm2Smooth);
phi2=V2-max(V2(:))/2;
%%
for i = 1:5
    phi2 = ac_ChanVese_model(max(phi2(:))*(V2), phi2, smooth_weight, image_weight, delta_t, 1); 
    
    if exist('h','var') && all(ishandle(h)),delete(h); end
  %  iso=isosurface(phi2);
 %  h = patch(iso,'facecolor','w');  axis equal;  view(3); 
  [C,h]=contour(phi2(:,:,midZ),[0,0],'w');axis equal;
  
  %  set(gcf,'name', sprintf('#iters = %d',i));
 %   drawnow; 
end
%%

phiBW=(phi2>0);
% for i=1:size(phi2,3);
%     phiBW(:,:,i)=bwconvhull(phiBW(:,:,i),'objects');
% end

[~,midPlane]=max(squeeze(sum(sum(double(phi2>0),1),2)));
phiBWMid=phiBW(:,:,midPlane);

%%
%phiBWMidBoundary= bwboundaries(phiBWMid);
[phiBWMidBoundary(:,1),phiBWMidBoundary(:,2)]= find(phiBWMid);
%%
%phiBWMidBoundary=cell2mat(phiBWMidBoundary);
shp = alphaShape(phiBWMidBoundary(:,1),phiBWMidBoundary(:,2));
aAlpha=40;
shp.Alpha=aAlpha;
bf = boundaryFacets(shp);
bf=[bf;bf(1,:)];
bfx=shp.Points(bf(:,1),1);
bfy=shp.Points(bf(:,1),2);


s=[0,cumsum(sqrt(diff(bfx).^2+diff(bfy).^2))']';
 L=max(s)*3;
    bfx=interp1(s,bfx,0:max(s)/L:max(s),'spline')';
    bfy=interp1(s,bfy,0:max(s)/L:max(s),'spline')';
bfxr=round(bfx);
bfyr=round(bfy);

phiBWMid(sub2ind(size(phiBWMid),bfxr,bfyr))=true;
phiBWMid=imfill(phiBWMid,'holes');

%%

BW3 = bwmorph(phiBWMid,'thin',Inf);

%%

for iFrame = 2:numFrames
    lastFrame = read(vidObj, iFrame);
    lastFrame=lastFrame(rect(1):rect(1)+rect(3),rect(2):rect(2)+rect(4));

    lastFrame=normalizeRange(sum(double(lastFrame),3));

lastFrame=bpass(lastFrame,.5,[30,30]);

    
    imshow(lastFrame);
    hold on
    for iIterations=1:10
    phi2 = ac_ChanVese_model(max(phi2(:))*lastFrame, phi2, smooth_weight, image_weight, delta_t, 1); 
    
    if exist('h','var') && all(ishandle(h)),delete(h); end
 %  h = patch(iso,'facecolor','w');  axis equal;  view(3); 
  [C,h]=contour(phi2,[0,0],'w');axis equal;
  
    set(gcf,'name', sprintf('#iters = %d , #frame = %d',iIterations,iFrame));
    drawnow; 
    end
    contours(iFrame).C=C;
    hold off
end






%%
figure;
slice = [10,15,20,25,30,35,40,45];
for i = 1:8
    subplot(2,4,i); imshow(V(:,:,slice(i)),[]); hold on; 
    c = contours(phi(:,:,slice(i)),[0,0]);
    zy_plot_contours(c,'linewidth',2);
end