
% run for immobilized worms to create a static centerline. First use
% wormCL_viewer to click a good centerline for the first frame. Then run
% this program on the Behavior folder to copy the manual centerline for all
% frames.
display('Select a Brainscanner folder.')
try
behaviorFolder=uipickfiles('FilterSpec',mostRecent);
catch
    behaviorFolder=uipickfiles();
end
behaviorFolder = behaviorFolder{1}
behaviorFolder = [behaviorFolder filesep BehaviorAnalysis]
load([behaviorFolder filesep 'centerline'])
% make a copy of previous CL
save([behaviorFolder filesep 'CL_copy'] ,'centerline','eigenProj'...
    ,'wormcentered');
% copy first centerline entry --we assume this is the one you want
cltmp = centerline(:,:,1);
length = size(centerline)
centerline2 = repmat(cltmp, [1,1,length(end)]);
eigenWormFile='eigenWorms.mat';
load(eigenWormFile);
wormcentered=FindWormCentered(centerline2);
eigbasis=imresize(eigbasis,[size(eigbasis,1),size(wormcentered,1)]);
eigenProj=eigbasis*wormcentered;
%save new centerlines
display(['Saving data in ' behaviorFolder])
save([behaviorFolder filesep 'centerline'] ,'centerline2','eigenProj'...
,'wormcentered');
