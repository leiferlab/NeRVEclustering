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
endop
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
%% testing single 2d
plane



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
C=C';
C(C(:,1)==0,:)=[];
contours(iSlice).C=C;
end

%%
C1=[contours(1).C];
C2=[contours(10).C];

[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(C1, C2, 3, [4, 0.2, 0.01], [0.0000008, 0.0000008, 0.0000008], [0 0],1,1);

%%
plot([Transformed_M(:,1) C1(:,1)]',[Transformed_M(:,2) C1(:,2)]','r')
hold on
scatter(C1(:,1),C1(:,2),'b')
interpF.method='nearest';interpF.radius=50;
[imgw, imgwr, map] = tpswarp(V1, size(V1'), Transformed_M(1:20:end,:), C1(1:20:end,:), interpF);
%%

[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(C1, C2, 3, [4, 0.2, 0.01], [0.0000008, 0.0000008, 0.0000008], [0 0 0],1,1);



%% testing single3d plane

margin = 10; 
% phi = zeros(size(V)); 
% phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1; 
% phi = ac_reinit(phi-.5); 
% 

V2=normalizeRange(V);
phi=V2-max(V2(:))/2;
%%
phi2=(phi);
for i = 1:1
    phi2 = ac_ChanVese_model(max(phi2(:))*(V2), phi2, smooth_weight, image_weight, delta_t, 5); 
    
    if exist('h','var') && all(ishandle(h)),delete(h); end
  %  iso=isosurface(phi2);
 %  h = patch(iso,'facecolor','w');  axis equal;  view(3); 
 % [C,h]=contour(phi2,[0,0],'w');axis equal;
  
    set(gcf,'name', sprintf('#iters = %d',i));
    drawnow; 
end
contours(1).C=C;
hold off


%%



%%
figure;
slice = [10,15,20,25,30,35,40,45];
for i = 1:8
    subplot(2,4,i); imshow(V(:,:,slice(i)),[]); hold on; 
    c = contours(phi(:,:,slice(i)),[0,0]);
    h=zy_plot_contours(c,'linewidth',2);
end