function clusterWormStraigthening(dataFolder,nStart,nRange)
clusterFolder=['/scratch/tmp/jnguyen/' dataFolder];
load([clusterFolder filesep 'startWorkspace.mat']);

for iStack=nStart:(nStart+nRange-1)
   try 
    WormCLStraighten_5(clusterFolder,destination,vidInfo,...
    alignments,[],Vtemplate,vRegion,zOffset,iStack,side,0);

   catch me
       me
   end
end
