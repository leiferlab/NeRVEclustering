function imageAll=zProjectDatFile(bgFile,imsize,method,skip)

%loads dat file and calculates the average intensity of each pixel. method
%is not yet available. image can be any size. possible methods are 'mean', 
%'sum','square', 'min', and 'max'
if nargin==0
bgFile=uipickfiles;
end
if isempty(bgFile);
    bgFile=uipickfiles;
end
if iscell(bgFile)
    bgFile=bgFile{1};
end

if nargin<2
    [row, col]=getdatdimensions(bgFile);
    imsize=[row col];
end
 if nargin<3
     method='mean';
 end
if nargin<4
    skip=1;
end


Fid=fopen(bgFile);

status=fseek(Fid,0,1);
lastFrame=ftell(Fid)/(prod(imsize)*2);
status=fseek(Fid,0,-1);
progressbar(0);
if  strcmp(method,'min')
    imageAll=Inf;
else
    imageAll=0;
end
for i=1:skip:lastFrame
progressbar(i/lastFrame);
  pixelValues=fread(Fid,prod(imsize),'uint16',0,'l');

squareImage=(reshape(pixelValues,imsize(1),imsize(2)));
if strcmp(method,'sum') || strcmp(method,'mean')
    imageAll=imageAll+double(squareImage);
    elseif strcmp(method,'square')
        imageAll=imageAll+double(squareImage).^2;
    elseif strcmp(method,'max');
    imageAll=max(imageAll,double(squareImage));
    elseif strcmp(method,'min')
    imageAll=min(imageAll,double(squareImage));
end

end
if strcmp(method,'mean') ||strcmp(method,'square')
imageAll=imageAll/lastFrame*skip;
end
