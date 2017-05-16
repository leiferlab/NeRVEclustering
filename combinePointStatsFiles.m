function pointStats=combinePointStatsFiles(psFolder,psLength)
% creates a complete pointStats variable from an analysis folder which is
% created by the WormStraightening programs. Adjusted to work on the
% cluster, DOES NOT OVERWRITE EXISTING POINTSTATS FILES

% output is a saved pointStats in the the dataFolder,
%       pointStats is an array of structures, with the ith structure
%       coresponding to the ith volume in the recording. 

if nargin==0
    dataFolder=uipickfiles;
    dataFolder=dataFolder{1};
end
if nargin==1
    dataFolder=psFolder;
end

psFolder=dir([dataFolder filesep 'CLstraight*']);
psFolder=[dataFolder filesep psFolder(end).name];
pList=dir([psFolder filesep 'pointStats*']);
if nargin==1
    psLength=length(pList);
end

pointStats=repmat(struct(),1,psLength);
%progressbar(0);
for iFile=1:length(pList)
    %   progressbar(iFile/length(pList));
    if ~mod(iFile,100)
        display([ 'Completed frame ' num2str(iFile) ' of ' num2str(psLength)])
    end
    
    idx=str2double(pList(iFile).name(11:15));
    %addint trycatch loop because matfile occasionally is corrupted during
    %saving
    try
        input=load([psFolder filesep pList(iFile).name]);
        input=input.pointStats;
        good_vol_flag= length(input.straightPoints)<200;
    catch me
        fprintf([me.message '\n'])
        good_vol_flag=0;
    end
    
    if good_vol_flag
        pointStats(idx).stackIdx=input.stackIdx;
        pointStats(idx).straightPoints=input.straightPoints;
        pointStats(idx).rawPoints=input.rawPoints;
        pointStats(idx).pointIdx=input.pointIdx;
        pointStats(idx).Rintensities=input.Rintensities;
        pointStats(idx).Volume=input.Volume;
    else
        pointStats(idx).stackIdx=[];
        pointStats(idx).straightPoints=[];
        pointStats(idx).rawPoints=[];
        pointStats(idx).pointIdx=[];
        pointStats(idx).Rintensities=[];
        pointStats(idx).Volume=[];
    end
end

%save final pointStats File
save([dataFolder filesep 'PointsStats'],'pointStats');

%write to status file
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
if contains(hostname,'della')
    Fid=fopen([dataFolder filesep 'status.txt'],'a');
    status=[datestr(datetime('now')) ':Finished straightening \n'];
    fprintf(Fid,status);
    fclose(Fid);
end
