function imOut=pixelIntensityCorrection(imIn)
% function loads pixelcalibration mat file and uses it to correct the
% intensities of imIn on a pixel by pixel basis. On the first call, it
% loads the calibration data and saves it to the window to save time on
% future calls. 


pixM=getappdata(0,'pixM');
pixB=getappdata(0,'pixB');
pixM(pixM<.7)=.7;
pixM(pixM>1.5)=1.5;
if isempty(pixM) || isempty(pixB)
try
    %try to load pixel calibration file, if not present, user picks a file
    load('Y:\CommunalCode\3dbrain\PixelCalibration\');
catch
calibration=uipickfiles;

load(calibration{1});
end
    setappdata(0,'pixM',pixM);
    setappdata(0,'pixB',pixB);
end

if isempty(pixM) || isempty(pixB)
    error('pixel stats not found');
end

%pixM(abs(pixM-1)>.2)=1;
imOut=(imIn-pixB).*pixM;
imOut(imOut<0)=0;
