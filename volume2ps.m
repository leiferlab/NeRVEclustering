function pointStats=volume2ps(varargin)
% volume2ps takes a volume image of the worm brain, usually already
% straightened. It runs the segmentation algorithm based on the eigenvalues
% of the hessian matrix, and then calculates shape metricds like brightness
% and volume for each object, storing the results in a pointStats structure
% commonly used throughout the wholebrain imaging pipeline. THERE IS A MAX
% NUMBER OF OBJECTS DETECTED SET AT 200. If this is exceeded, the program
% will return empty matricies, as too many points is a signal that the
% image is bad. 

% JN 20170515

%Inputs
% pointStats=volume2ps(V)
%       V - an input 3D fluorecent image to be segmented.This will use
%       default segmentation options and will run its own round of
%       smoothing
%
% pointStats=Volume2ps(V,options)
%       options - stucture of segmentation options to be used. The user can
%       provide as many or as few options as they wish. Default values will
%       be used when no field is provided.
%
% pointStats=Volume2ps(V,options,Vsmooth)
%       Vsmooth - optional input for a pre-smoothed version of the volume.
%       If this is not provided, the code will run its own bandpass filter
%       on V. Providing Vsmooth makes the code run faster. 




V=varargin{1};

%default options for segmentation
options.method='invdist';
options.radius=20;
options.power=1;
options.thresh1=.05;
options.minObjSize=50;
options.maxObjSize=400;
options.minSphericity=.80;
options.filterSize=[10 10 4];
options.power=1;
options.prefilter=0;
options.hthresh=0;

%load user provided inputs inf any
if nargin>=2
    options_provided=fieldnames(varargin{2});
    for i=1:length(options_provided)
        options.(options_provided{i})=varargin{2}.(options_provided{i});
    end
end


%run segmentation
if nargin==3
    Vsmooth=varargin{3};
    [wormBW2,~]=WormSegmentHessian3dStraighten(V,options,Vsmooth);
    
else
    [wormBW2,~]=WormSegmentHessian3dStraighten(V,options);
    
end

% option, to use presmoothed version, much faster but may not be a s good
% do segmentation on straightened worm

imsize=size(V);

%% Filter out noise at top and bottom of stack
% bpass filters tend to add noise at the top and the bottom of the
% stack, we want to remove this. This will be done by projecting along
% the x and y of the image mask to get a plot of number of pixels vs Z.
% We expect to find a peak in the number of pixels close to the ends of
% the stack.

%project along x and y
BWplot=(squeeze(sum(sum(wormBW2,1),2)));
%smooth
BWplot=smooth(BWplot,20);
%find peaks on  the projection plot
[~,peak_loc]=findpeaks(BWplot);
%find lowest and highest peaks, corresponding to the peaks in noise
%near the edges
endpts=peak_loc([1,length(peak_loc)]);
%find the trough by finding peaks of -BWplot,
[~,locs]=findpeaks(-BWplot);
%find the ones closes the top and bottom of the worm, this will define
%deltion zones between 1 and botpoint1, and botpoints2 and the end of
%the stack
botpoint1=locs((locs>endpts(1)));
botpoint2=locs((locs<endpts(2)));
% if no troughs are found, great, just use the top and bottom of the
% stack (this will end up doing nothing in the end)
if isempty(botpoint1);botpoint1=1;end
botpoint1=botpoint1(1);

if isempty(botpoint2);botpoint2=imsize(3);end
botpoint2=botpoint2(end);
% if the troughs founds are more than a quarter away from the edges,
% then there probably isnt any noise artifact at the edges, so move
botpoint1(botpoint1>imsize(3)*1/4)=1;
botpoint2(botpoint2<imsize(3)*3/4)=imsize(3);


cc=bwconncomp(wormBW2,6);

%kill of objects that touch the deletion zone
%function that take linear pixel idx and returns Z
zindexer=@(x,s) x./(s)+1;
%apply zindexer to each of the objects from bwconncomp
objectZ=cellfun(@(x) zindexer(x,imsize(1)*imsize(2)),cc.PixelIdxList...
    ,'uniformoutput',0);

%label objects as bad if theh touch the deletion zones
badRegions_bot=cellfun(@(x) any(x<=botpoint1),objectZ);
badRegions_top=cellfun(@(x) any(x>=botpoint2),objectZ);
badRegions=(badRegions_bot | badRegions_top)';

% remove bad objects by setting all their pixels to false
wormBW2(cell2mat(cc.PixelIdxList(badRegions)'))=false;
cc.PixelIdxList=cc.PixelIdxList(~badRegions);
cc.NumObjects=nnz(~badRegions);


%% compile results into pointStats structure for saving
%hard cap at 200 neurons, occasionally  you have big fails that produce
%many hundreds of points. this is bad, just blank everything if this
%happens
cc=bwconncomp(wormBW2,6);

%get results from binary mask using region props
stats=regionprops(cc,V,'Centroid','MeanIntensity',...
    'Area');

%if there are too many points, get rid of some of the dimmer ones
if cc.NumObjects>150
    intensities=[stats.MeanIntensity];
    intensities_sort=sort(intensities,'descend');
    badRegions=intensities<intensities_sort(150);
    
    
    % remove bad objects by setting all their pixels to false
    wormBW2(cell2mat(cc.PixelIdxList(badRegions)'))=false;
    cc.PixelIdxList=cc.PixelIdxList(~badRegions);
    cc.NumObjects=nnz(~badRegions);
    stats=stats(~badRegions);
end

%get roi mean intensities
intensities=[stats.MeanIntensity]';
%get roi centroids
P=cell2mat({stats.Centroid}');

%flip x and y
P(:,[1 2])=P(:,[2 1]);

%get roi volumes
Areas=[stats.Area]';
%centroids from unstraightened coordinate system using
%interpolation.


% load up all of these into pointstats
pointStats.straightPoints=P(:,1:3);
pointStats.pointIdx=(1:cc.NumObjects)';
pointStats.Rintensities=intensities;
pointStats.Volume=Areas;

% save the binary mask

pointStats.baseImg=logical(wormBW2);

