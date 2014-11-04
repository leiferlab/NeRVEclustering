%select file name

movieFile=uipickfiles('FilterSpec','Y:\PanNeuronal\');

%% load vidObj using VideoReader
vidObj = VideoReader(movieFile{1});
lastFrame = read(vidObj, inf);
 numFrames= vidObj.NumberOfFrames;

%%
smooth_weight = .001; 
image_weight = 100; 
delta_t = 4; 

%%
firstFrame=read(vidObj, 1);
firstFrame=im2double(firstFrame(:,:,1));
firstFrame=normalizeRange(firstFrame);
phi=firstFrame-.5;
%%
contours.C=[];
contours=repmat(contours,1,numFrames);
nIterations=20;
m=max(phi(:));
%%
parfor iFrame = 1:numFrames
    tic
    try
    lastFrame = read(vidObj, iFrame);
lastFrame=im2double(lastFrame(:,:,1));
lastFrame=normalizeRange(lastFrame);
phi=lastFrame-.5;
m=max(phi(:));
%lastFrame=bpass(lastFrame,.5,[30,30]);
    phi = ac_ChanVese_model(m*lastFrame, phi, smooth_weight, image_weight, delta_t, nIterations); 
    C=contourc(phi,[0 0]);
    contours(iFrame).C=C';
    display(['Finished frame ' num2str(iFrame) ' in ' num2str(toc) ' s']);
    catch me
        contours(iFrame).C=[];
        display(['Error frame ' num2str(iFrame) ' in ' num2str(toc) ' s']);

    end
    
 end


% match contours with neighbours
Cin2=contours(2:numFrames);
Ctransform=repmat(struct,1,numFrames);
groupSize=100;
for iGroup=1:numFrames/groupSize;
    group=(1+(iGroup-1)*groupSize:iGroup*groupSize);
    group=group(group>0 |group<numFrames);
parfor i=group
    try
C1=contours(i).C;
C2=Cin2(i).C;

C1=C1(C1(:,1)~=0,:);
C2=C2(C2(:,1)~=0,:);
tic
[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(C1, C2, 2, [.1, 0.01], [0.000008, 0.00008, 0.0000008], [0 0],1,0);
Ctransform(i).Transformed_M=Transformed_M;
Ctransform(i).multilevel_ctrl_pts= multilevel_ctrl_pts;
Ctransform(i).multilevel_param=multilevel_param;

display(['Finished ' num2str(i) ' fit in '  num2str(toc) ' s' ])
    catch
        display(['Error ' num2str(i) ])
    end
end

end
