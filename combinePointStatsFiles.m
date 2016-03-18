function pointStats=combinePointStatsFiles(psFolder,psLength)
% creates a complete pointStats variable from an analysis folder which is
% created by the WormStraightening programs.
if isempty(psFolder)
    psFolder=uipickfiles;
    psFolder=psFolder{1};
end
pointStats=repmat(struct(),1,psLength);
pList=dir([psFolder filesep 'pointStats*']);

progressbar(0);
for iFile=1:length(pList);
    progressbar(iFile/length(pList));
    idx=str2double(pList(iFile).name(11:15));
    input=load([psFolder filesep pList(iFile).name]);
    input=input.pointStats;
    pointStats(idx).straightPoints=input.straightPoints;
    pointStats(idx).rawPoints=input.rawPoints;
    pointStats(idx).stackIdx=input.stackIdx;
    pointStats(idx).pointIdx=input.pointIdx;
    pointStats(idx).Rintensities=input.Rintensities;
    pointStats(idx).Volume=input.Volume;
end