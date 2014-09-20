function imageAll=zProjectDatFile(bgFile,method)

%loads dat file and calculates the average intensity of each pixel. method
%is not yet available. image must be 1024x1024, and be uint16. 
% WARNING: I stick a transpose after the reshape.
if nargin==0
bgFile=uipickfiles;
bgFile=bgFile{1};
end


Fid=fopen(bgFile);

status=fseek(Fid,0,1);
lastFrame=ftell(Fid)/(1024^2*2);
status=fseek(Fid,0,-1);
progressbar(0);
for i=1:lastFrame
progressbar(i/lastFrame);
  pixelValues=fread(Fid,1024^2,'uint16',0,'l');

squareImage=(reshape(pixelValues,1024,1024)');
if i==1
    imageAll=double(squareImage);
else
    imageAll=imageAll+double(squareImage);
end

end
imageAll=imageAll/lastFrame;