function WormEigenBehavior(dataFolder)
% WormEigenBehavior takes an input as a worm data folder,the user selects a
% centerline file in that folder and the program creates the behavior
% analysis folder with the centerline , eigenprojections, and wormcentered.
% Data can then be loaded with the loadCLBehavior function. 
CLfile=uipickfiles('filterspec',dataFolder);
centerlineFile=load(CLfile{1});
fieldNames=fieldnames(centerlineFile);
clIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'line'));

    centerline=centerlineFile.(fieldNames{clIdx(1)});
%create wormcentered coordinate system
wormcentered=FindWormCentered(centerline);
%project onto eigen basis
load('Y:\CommunalCode\3dbrain\eigenWorms.mat')
if size(eigbasis,2)~=size(wormcentered,1)
    eigbasis=imresize(eigbasis,[size(eigbasis,1),size(wormcentered,1)]);
    
end
eigenProj=eigbasis*wormcentered;

%save outputs into behaviorAnalysis folder
behaviorFolder=[dataFolder filesep 'BehaviorAnalysis'];
mkdir(behaviorFolder);
offset=0;
save([behaviorFolder filesep 'centerline'] ,'centerline','eigenProj'...
    ,'wormcentered','offset');

