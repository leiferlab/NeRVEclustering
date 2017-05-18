% This is a tutorial script shows how to create neuron registration vectors
% for a single volume with a set of reference files.  

%% set up the path
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );

ver=version('-release');
ver=str2double(ver(1:4));
if ver<2017
    error('Please use a version of matlab newer than 2017');
end
    

if contains(hostname,'tigressdata')
    cd /tigress/LEIFER/communalCode/3dbrain/
    path(pathdef)
    fprintf('Adding files to path')
else
    disp(['This code is designed to work on tigressdata.'...
        ' You''re not currently on tigressdata so make sure you have the',...
        ' 3dbrain repo in your path!'])
end




%% select the folder Brainscanner folder, this will have all the things we need
% the folder has all the data from the analysis pipeline saved in various
% files. For the time being, all the file names are set so no additional
% inputs are required. 

dataFolder=uipickfiles('Prompt', 'Select the Brain folder', ...
    'FilterSpec','/tigress/LEIFER/PanNeuronal/testing_sets');
dataFolder=dataFolder{1};

%% First we'll load up the pointStats file within the dataFolder

ps_filename=[dataFolder filesep 'pointStatsNew.mat'];
pointStats=load(ps_filename);
pointStats=pointStats.pointStatsNew;

%there are several pointStats files in the folder, as they are modified
%along the analysis pipeline. This is the final one. Again, names should be
%changed when convenient.

%% the timing can also be found in heatData
heatFile=[dataFolder filesep 'heatData.mat'];
heatData=load(heatFile);

%the timing for each volume is present in heatData.hasPointsTime. 

%% notes on timing and relation to other structures

%from the hiResImages tutorial
[~,~,hiResData]=tripleFlashAlign(dataFolder);

%recall hiResData.stackIdx shows all the frames that belong to a particular
%volume. The number in stackIdx matches the array position in pointStats,
%so the frames with hiResData.stackIdx==100 are used to create
%pointStats(100). 


%% looking ah the fields

% pointStats has information from each volume, derived from the
% straightened images.This comes from both segmentation and error
% correction.

%to get data from the 100th volume, just pull the 100th structure in the
%array. 
target_volume=199;

pointStats(target_volume)


%% Straightened coordinates


%plotting straightened cooridinates from one volume. The coordinates are in
%units of himag pixels, with each voxel being a cube (no rescaling required
%between xy and z. 

ps=pointStats(target_volume);

scatter3(ps.straightPoints(:,1),...
    ps.straightPoints(:,2),...
    ps.straightPoints(:,3))
title(['Example all straightened neuron coordinates, Vol: '...
    num2str(target_volume)])
axis equal


% the points displayed here are ordered in how they were found during
% segmentation, so row 1 in volume 1 does not correspond to row 1 in volume
% 2. You will see that straightPoints and trackIdx will be longer than the
% other fields in the structure. That is because these fields include
% contributions from error correction while the other ones do not.

%% looking at tracked points
%not all of the coordinates in straightPoints represent neurons that are
%tracked throughout the recording. Only the ones that have a non-nan value
%of trackIdx will have labels across frames. The value of trackIdx is what
%identifies the neurons across frames,

target_volume=199;
ps=pointStats(target_volume);

%get the points that are tracked
present_idx=~isnan(ps.trackIdx);
present_points=ps.straightPoints(present_idx,:);
%get the index of the present points
track_idx=ps.trackIdx(present_idx);


%reorganize the points in the proper order
tracked_points=zeros(size(present_points));
tracked_points(track_idx,:)=present_points;

% now the coordinates in tracked_points are ordered, so after following
% this procedure for each volume, neuron one will be the first row for
% every volume. 


scatter3(tracked_points(:,1),...
    tracked_points(:,2),...
    tracked_points(:,3))
xlabel('x');ylabel('y');zlabel('z')
title(['Example tracked neuron coordinates from volume' num2str(target_volume)])
axis equal


%% lets do it again with another volume
% to visualize the matching between volumes. 

target_volume2=299;
ps=pointStats(target_volume2);

%get the points that are tracked
present_idx=~isnan(ps.trackIdx);
present_points=ps.straightPoints(present_idx,:);
%get the index of the present points
track_idx=ps.trackIdx(present_idx);


%reorganize the points in the proper order
tracked_points2=zeros(size(present_points));
tracked_points2(track_idx,:)=present_points;

%scatter points from the first volume
scatter3(tracked_points(:,1),...
    tracked_points(:,2),...
    tracked_points(:,3),'xr')


hold on
%scatter points from the other volume
scatter3(tracked_points2(:,1),...
    tracked_points2(:,2),...
    tracked_points2(:,3),'ob')

%plot lines between corresponding neurons
plot3([tracked_points(:,1),tracked_points2(:,1)]',...
    [tracked_points(:,2),tracked_points2(:,2)]',...
    [tracked_points(:,3),tracked_points2(:,3)]','g')
title('Neuron matching in different frames in straightened coordinate system')
axis equal
hold off
legend(['Vol: ' num2str(target_volume)],['Vol: ' num2str(target_volume2)]);

% in this order, the row 1 will also correspond witht he first row of all
% of the signal matrices in heatData. tracked_points should have the same
% number of rows as the signal fields in heatdata;


%% other fields

%there are various other fields in each pointStats structure that were
%required for tracking. Rintensities and Volume were used for non-rigid
%pointset registration to represent the coordinates as a Gaussian mixture.
%Raw points are the coordinate locations in the original frame of
%reference, but doesnt include points from error correction so this isn't
%very useful at the moment. trackMatrix comes is the registration vector
%for each of the neurons found before error correction. trackWeights is the
%value of the projection of each registration vector onto the average
%vector for the cluster that neuron was assigned to.


