function [xCent,yCent,majAxis,minAxis,orientation] = xyCentSize(subImageIn);
% given a 'small' image of a bright object on a zero based background,
% find the the center of mass in x and y
% find the major/minor axes and orientation of an elliptical contour at 50% of max height
% 

%% center of mass
totalIntensity = sum(subImageIn(:));
rowSum = sum(subImageIn,2);
colSum = sum(subImageIn,1);

xCent = sum((1:size(subImageIn,2)).*colSum)./totalIntensity;
yCent = sum((1:size(subImageIn,1))'.*rowSum)./totalIntensity;
%%
% % [xx,yy] = meshgrid(1:7,1:7);
% % [xi,yi] = meshgrid(1:0.6:7,1:7);
% z2 = -3:1:3;
% kernelZ = 1;
% 
% kernel = exp(-(z2.^2)/(2*(kernelZ^2)));
% kernel = kernel./sum(kernel(:));
% size1 = size(kernel,2);
% kernel2 = repmat(kernel,[size1,1]);
% kernel3 = kernel2.*kernel2';
% 
% [xx,yy] = meshgrid(1:size1,1:size1);
% [xi,yi] = meshgrid(1:0.6:size1,1:size1);
% subImageIn = interp2(xx,yy,kernel3,xi,yi);
% % subImageIn = kernel3;
% 
% subImageIn = imrotate(subImageIn,35,'Bicubic');
% imshow(subImageIn,[],'InitialMagnification','fit');
% figure(gcf);
%% interpolated image
contourLevel = 0.5*max(subImageIn(:));
nInterp = 2;
subImageTest = interp2(subImageIn,nInterp,'Bicubic');
subImageTest = subImageTest>contourLevel;
propsOut = regionprops(single(subImageTest),'Centroid','Orientation','MajorAxisLength','MinorAxisLength');
orientation = propsOut.Orientation;
majAxis = propsOut.MajorAxisLength/(2^nInterp);
minAxis = propsOut.MinorAxisLength/(2^nInterp);



