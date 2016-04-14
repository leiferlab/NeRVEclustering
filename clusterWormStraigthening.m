function clusterWormStraigthening(dataFolder,nStart,nRange)
% calls worm straightening code if a startWorkspace.mat file is already
% created in the dataFolder being analyzed, program runs to straighten
% stacks nStart:nStart+nRange-1
straightenData=load([dataFolder filesep 'startWorkspace.mat']);

destination=straightenData.destination;
Vtemplate=straightenData.Vtemplate;
zOffset=straightenData.zOffset;
side=straightenData.side;
lastOffset=straightenData.lastOffset;
vidInfo=straightenData.vidInfo;
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
ctrlPoints=[];
display(dataFolder)
imageFolder2=[dataFolder filesep destination];

for iStack=nStart:(nStart+nRange-1)
    fileName2=[imageFolder2 filesep 'image' num2str(iStack,'%3.5d') '.tif'];
    fileName3=[imageFolder2 filesep 'pointStats' num2str(iStack,'%3.5d')];
    if ~exist(fileName2) && ~exist(fileName3)
        WormCLStraighten_11(dataFolder,destination,vidInfo,...
            alignments,ctrlPoints,Vtemplate,zOffset,iStack,side,lastOffset,0);
    end
end
