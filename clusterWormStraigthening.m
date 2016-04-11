function clusterWormStraigthening(dataFolder,nStart,nRange)
% calls worm straightening code if a startWorkspace.mat file is already
% created in the dataFolder being analyzed, program runs to straighten
% stacks nStart:nStart+nRange-1
load([dataFolder filesep 'startWorkspace.mat']);
display(dataFolder)

for iStack=nStart:(nStart+nRange-1)
   try 
       WormCLStraighten_11(dataFolder,destination,vidInfo,...
            alignments,ctrlPoints,Vtemplate,zOffset,iStack,side,lastOffset,show);
   catch me
       for i=1:length(me.stack)
       me.stack(i)
       end
       save([dataFolder filesep 'Error' num2str(iStack)],'me')
   end
end
