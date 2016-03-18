function [fiducialPoints,z2ImageIdxOffset]=loadFiducialPoints(dataFolder)
if nargin>=1
%% load Fiducials file
fiducialFileName=dir([dataFolder filesep '*iducial*']);
fiducialFileName={fiducialFileName.name}';
else
    fiducialFileName=[];
end

if length(fiducialFileName)~=1
        display('Select model file');

    fiducialFileName=uipickfiles('FilterSpec',dataFolder);
    fiducialFile=load(fiducialFileName{1});
    fiducialPoints=fiducialFile.fiducialPoints;
else
    try
    fiducialFile=load([dataFolder filesep fiducialFileName{1}]);
    catch
    display('Select model file');
    fiducialFileName=uipickfiles('FilterSpec',dataFolder);
    fiducialFile=load(fiducialFileName{1});
    end
    
    fiducialPoints=fiducialFile.fiducialPoints;

end
%%
s1=max(cellfun(@(x) size(x,1),fiducialPoints));
s2=max(cellfun(@(x) size(x,2),fiducialPoints));
for i=1:length(fiducialPoints);
    fiducialPointsTemp=fiducialPoints{i};
    stemp=size(fiducialPointsTemp);
    if stemp(2)-s2<0
        fiducialPointsTemp=[fiducialPointsTemp cell(stemp(1),s2-stemp(2))];
    end
    if stemp(1)-s1<0
        fiducialPointsTemp=[fiducialPointsTemp; cell(s1-stemp(1),s2)];
        
    end
    fiducialPoints{i}=fiducialPointsTemp;
    
    
    
end
if nargout==2
%
try
    z2ImageIdxOffset=fiducialFile.timeOffset;
catch
   timeFile= dir([fileparts(fiducialFileName{1}) filesep '*ffset*']);
   if isempty(timeFile)
    timeFile=uipickfiles('FilterSpec',fileparts(fiducialFileName{1}));
    timeFile=timeFile{1};
   else
       timeFile=[fileparts(fiducialFileName{1}) filesep timeFile.name];
   end
   
    timeFile=load(timeFile);
    z2ImageIdxOffset=timeFile.timeOffset;

end
end

