function clusterWormStraigthening(dataFolder,nStart,nRange)
[~,dataFolder]=fileparts(dataFolder);

clusterFolder=['/scratch/tmp/jnguyen/' dataFolder];
display(clusterFolder)
load([clusterFolder filesep 'startWorkspace.mat']);
display(clusterFolder)

for iStack=nStart:(nStart+nRange-1)
   try 
    WormCLStraighten_5(clusterFolder,destination,vidInfo,...
    alignments,[],Vtemplate,vRegion,zOffset,iStack,side,0);

   catch me
       for i=1:length(me.stack)
       me.stack(i)
       end
       save([clusterFolder filesep 'Error' num2str(iStack)],'me')
   end
end
