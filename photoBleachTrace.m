%% photoBleachTrace for whole field of view


display('Select datFile');
datFile=uipickfiles('FilterSpec','O:\');
datFile=datFile{1};
dataFolder=fileparts(datFile);
timeStep=.005; % 5 ms between frames
imFlash=datTimeTrace(datFile,1200,600);
% fit result to exponential (linear fit to semi-log)
p=polyfit((50000:length(imFlash))*timeStep,log(imFlash(50000:end)),1);
timeConstant=-1/p(1);
save([dataFolder filesep 'photoBleach'],'imFlash','timeConstant');
display(['Decay Constant is : ' num2str(timeConstant) 's'])