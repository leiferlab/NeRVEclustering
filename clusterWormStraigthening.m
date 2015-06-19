function clusterWormStraigthening(dataFolder,nStart,nRange)
clusterFolder=['/scratch/tmp/jnguyen/' dataFolder];
load([clusterFolder filesep 'startWorkspace.mat']);

for iStack=nStart:(nStart+nRange-1)
    
[~,pointStats,~,~]=...
    WormCLStraighten_5(clusterFolder,destination,vidInfo,...
    alignments,[],Vtemplate,vRegion,zOffset,iStack,side,0);

save([dataFolder filesep  'pointStats' num2str(iStack,'%3.5d')],'pointStats');

end
