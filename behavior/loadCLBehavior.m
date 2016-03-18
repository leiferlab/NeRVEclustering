function [centerline, offset,eigenProj, CLV,wormcentered]=loadCLBehavior(dataFolder)
% input is a folder dataFolder with: 
% "./BehaviorAnalysis/centerline.mat
% made to be a little bit flexible with file names

%loadCLBehavior takes in centerline data created by WormEigenBehavior.m
%into workspace. It also calculates the sign of the centerline velocity

%select centerline file
cldir=dir([dataFolder filesep 'BehaviorAnalysis' filesep 'centerlin*']);

centerlineFile=load([dataFolder filesep 'BehaviorAnalysis' filesep cldir(1).name]);
%centerline file should be a structure with cetnerline, wormcentered, and
%eigProj.
fieldNames=fieldnames(centerlineFile);
eigIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'eig'));
%create variables from loaded structures
for iField=1:length(fieldNames)
    currentStr=fieldNames{iField};
if strfind(currentStr,'wormcen')
    wormcentered=centerlineFile.(currentStr);
elseif strfind(currentStr,'line')
    centerline=centerlineFile.(currentStr);
end
end

if isfield(centerlineFile,'offset');
offset=centerlineFile.offset;
else
    offset=0;
end
% create other variables if they are called for in outputs

if nargout>2
    %if eigenProj is in the loaded structure, use that, otherwise, see if
    %the file is there.
if any(eigIdx)
    eigenProj=centerlineFile.(fieldNames{eigIdx});
else
    eigenProj=load([dataFolder filesep 'BehaviorAnalysis' filesep 'eigenProj']);
    fieldNames=fieldnames(eigenProj);
eigenProj=eigenProj.(fieldNames{1});
end

   %calculate sign of the worm velocity by looking at the two derivatives. 
wormcentered2=smooth2a(wormcentered,5,31);
[wormcenteredx,wormcenteredy]=gradient(wormcentered2);
CLV=sum((wormcenteredy.*wormcenteredx));

end