function [newImage,newBasisX,newBasisY]=wormStraightening(CL,imageIn,windowSize,pixSize)
%wormStraightening takes centerline data CL and an image and interpolates a
%curvilinear coordinate system around it.
if nargin==3
    pixSize=1;
end


windowSpan=-windowSize:pixSize:windowSize;
%find t vectors
CL2(:,1)=interp1(CL(:,1),1:pixSize:length(CL));
 CL2(:,2) =interp1(CL(:,2),1:pixSize:length(CL));
 CL=CL2;
 
 CL=[smooth(CL(:,1),20),smooth(CL(:,2),20)];
[~,perpSlopes]=gradient(CL,15/pixSize);
%create perpendicular n vectors
perpSlopes=normr([-perpSlopes(:,2),perpSlopes(:,1)]);


newBasisX=perpSlopes(:,1)*windowSpan;
newBasisY=perpSlopes(:,2)*windowSpan;

newBasisX=bsxfun(@plus,newBasisX,CL(:,1));
newBasisY=bsxfun(@plus,newBasisY,CL(:,2));

% imagesc(imageIn)
% hold on
% plot(CL(:,2),CL(:,1),'b');
% plot(newBasisY',newBasisX','black')

if ~isempty(imageIn)
newImage=interp2(sum(imageIn,3)',newBasisX,newBasisY);
else
    newImage=[];
end






