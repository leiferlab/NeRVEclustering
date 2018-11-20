function [centerline,eigenProj, CLV,wormcentered]=loadCLBehavior(dataFolder)
% input is a folder dataFolder with: 
% "./BehaviorAnalysis/centerline.mat
% made to be a little bit flexible with file names

%if any of the data is missing, calculate missing data and resave the file
    resaveFlag=0;


%loadCLBehavior takes in centerline data created by WormEigenBehavior.m
%into workspace. It also calculates the sign of the centerline velocity

%select centerline file
cldir=dir([dataFolder filesep 'BehaviorAnalysis' filesep 'centerlin*']);
if isempty(cldir)
    error('CL:missing', ...
   ['Centerline file not found! Ensure that BrainScanner has a'...
' "BehavioralAnalysis" folder and a "centerline.mat" file in it' ...
' or calculate centerlines']);
end

fileName=[dataFolder filesep 'BehaviorAnalysis' filesep cldir(1).name];
centerlineFile=load(fileName);
%centerline file should be a structure with cetnerline, wormcentered, and
%eigProj.
fieldNames=fieldnames(centerlineFile);
eigIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'eig'));
lineIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'line'));
cenIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'wormcen'));

centerline=centerlineFile.(fieldNames{lineIdx});

if any(cenIdx)
    wormcentered=centerlineFile.(fieldNames{cenIdx});
else
    wormcentered=FindWormCentered(centerline);
    resaveFlag=1;
end

if any(eigIdx)
    eigenProj=centerlineFile.(fieldNames{eigIdx});
else
    eigenWormFile='eigenWorms.mat';
    eigData=load(eigenWormFile);
    eigbasis=eigData.eigbasis;
    eigbasis=imresize(eigbasis,[size(eigbasis,1),size(wormcentered,1)]);
    eigenProj=eigbasis*wormcentered;
    resaveFlag=1;
end

   %calculate sign of the worm velocity by looking at the two derivatives. 
wormcentered2=smooth2a(wormcentered,5,31);
[wormcenteredx,wormcenteredy]=gradient(wormcentered2);
CLV=sum((wormcenteredy.*wormcenteredx));

if resaveFlag
    offset=0;
    display('Eigenprojections missing, calculating and resaving file');
    save(fileName,'centerline','wormcentered','eigenProj','offset');
end
