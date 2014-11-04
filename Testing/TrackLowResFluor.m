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
for iFrame=1:20:min(5000,numFrames)
    
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

%% find image gradients for tracking?


            

%%
cline_para.iterations=800;
cline_para.group_iterations=40;%40
cline_para.it_loop=80;%70
cline_para.pts=40; %pts in circle (needs to be divisible by 4)
cline_para.iterations3d=300;
if rem(cline_para.pts,4)~=0;
    cline_para.pts=cline_para.pts+4-rem(cline_para.pts,4);
end
cline_para.radius=20; %initial radius

cline_para.threshold=.1;
cline_para.stiff=0.05; %.0005
cline_para.stiff3d=1; %.1
cline_para.alpha =.005 ;
cline_para.beta = .5;
cline_para.alpha3d =00;
cline_para.beta3d =10;
cline_para.zalpha3d=0;
cline_para.zbeta3d=5;

cline_para.CLbeta=0.001;
cline_para.CLalpha=0.001;
cline_para.inter_frame_viscos_factor=0; %0
cline_para.kappa3d = 0; %1  %energy term
cline_para.kappa=100;
cline_para.gamma=40; %step size
cline_para.gradient_force =.5;%5
cline_para.interpolate_gradient_flag =0;
cline_para.show_flag=10; %1 to show end of every frame, 2 to show every compelted slice
cline_para.movie_flag=0;
cline_para.gradFlag=1;
cline_para.groupIterations=300;
cline_para.objsize=3;
cline_para.repulsion=20;

cline_para.loopSpringForce=1;
%%

lowResTest=smooth2a(lowResTest,1,1);

            %% parameters 
clear cline_initial2
clear cline_initial initialLoops

               h=imagesc(imStackSmooth(:,:,1));
            numLoops=inputdlg('How many loops?');
            numLoops=str2double(numLoops);
            %%
            for iLoop=1:numLoops

            cline_initial=imfreehand('closed', true);
            cline_initial=getPosition(cline_initial);
            cline_initial=[cline_initial;cline_initial(1,:)];
            s=[0;cumsum(sqrt(diff(cline_initial(:,1)).^2+diff(cline_initial(:,2)).^2))];
 
    cline_initial2(:,1)=interp1(s,cline_initial(:,1),0:1:max(s),'spline');
    cline_initial2(:,2)=interp1(s,cline_initial(:,2),0:1:max(s),'spline');
   cline_initial=cline_initial2; 
   cline_initial=cline_initial(1:end-1,:);
   
   initialLoops(iLoop).xyzs=cline_initial;
   initialLoops(iLoop).CM=mean(cline_initial);
   clear cline_initial2
            end
        
 cline_para.distance0=sqrt(sum(diff(cell2mat({initialLoops.CM}')).^2));

   
            
%%

cline_para.CLbeta=0.00;
cline_para.CLalpha=0.000001;
cline_para.kappa=10;

clineOutput = ActiveMultiLoopFit(imStackSmooth(:,:,:), cline_para,initialLoops); %gets center line and end pts

            
            

%%          


