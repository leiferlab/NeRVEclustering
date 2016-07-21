function makeHamSummary()

folderName=uigetdir;
tifNames=dir([folderName filesep '*.tif']);
progressbar(0)
for iImages=1:length(tifNames);
    [metaData(iImages),~]=getHamMetadata([folderName filesep tifNames(iImages).name]);
    progressbar(iImages/length(tifNames));
end
