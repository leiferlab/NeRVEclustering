function makePointStatsRef(dataFolder,nRef)
% create and save a pointStats mat file with only the volumes used as
% references. This saves the references used so that other volumes can
% match to the references consistently. This program can eventually be
% modified to select the references in a smarter way than just uniformly
% spaced. 

% Inputs:
%       filePath: the dataFolder containing the pointStats.mat file
%       nRef: the number of reference volumes to use

%Output:
%       a file pointStatsRef.mat is saved inside the dataFolder. 


PS_file=[dataFolder filesep 'pointStats.mat'];
load(PS_file)

%list of stacks presents
volList=find(cellfun(@(x) ~isempty(x),{pointStats.stackIdx}));
nVol=length(pointStats);

%the list of references to use
refList=unique(round(1:nRef/nVol:nVol));
refList=volList(refList);
PS_ref=pointStats(refList);
save([dataFolder filesep 'pointStatsRef.mat'],PS_ref)


