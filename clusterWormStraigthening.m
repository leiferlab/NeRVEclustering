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

for iStack=nStart:(nStart+nRange-1)
   try 
       WormCLStraighten_11(dataFolder,destination,vidInfo,...
            alignments,ctrlPoints,Vtemplate,zOffset,iStack,side,lastOffset,0);
   catch me
       for i=1:length(me.stack)
       me.stack(i)
       end
       save([dataFolder filesep 'Error' num2str(iStack)],'me')
   end
end
