function [newImage,newBasisX,newBasisY]=wormStraightening(CL,imageIn,windowSize,pixSize)
%wormStraightening takes centerline data CL and an image and interpolates a
%curvilinear coordinate system around it.
if nargin==3
    pixSize=1;
end

windowSpan=-windowSize:pixSize:windowSize;
%find t vectors
if ismatrix(CL)
    CL2(:,1)=interp1(CL(:,1,:),1:pixSize:size(CL,1))';
 CL2(:,2) =interp1(CL(:,2,:),1:pixSize:size(CL,1))';

else
CL2(:,1,:)=interp1(CL(:,1,:),1:pixSize:size(CL,1));
 CL2(:,2,:) =interp1(CL(:,2,:),1:pixSize:size(CL,1));
end
 CL=CL2;

 [x,dim,t]=size(CL);

 %CL=[smooth(CL(:,1),20),smooth(CL(:,2),20)];
[~,perpSlopes]=gradient(CL,15/pixSize);
%create perpendicular n vectors
perpSlopes=([-perpSlopes(:,2,:),perpSlopes(:,1,:)]);
perpSlopes=bsxfun(@rdivide,perpSlopes,sqrt(sum(perpSlopes.^2,2)));
perpX=perpSlopes(:,1,:);
perpY=perpSlopes(:,2,:);

newBasisX=reshape(perpX(:)*windowSpan,x,t,length(windowSpan));
newBasisX=permute(newBasisX,[1,3,2]);
newBasisY=reshape(perpY(:)*windowSpan,x,t,length(windowSpan));
newBasisY=permute(newBasisY,[1,3,2]);


newBasisX=bsxfun(@plus,newBasisX,CL(:,1,:));
newBasisY=bsxfun(@plus,newBasisY,CL(:,2,:));

% imagesc(imageIn)
% hold on
% plot(CL(:,2),CL(:,1),'b');
% plot(newBasisY',newBasisX','black')

if ~isempty(imageIn)
newImage=interp2(sum(imageIn,3)',newBasisX,newBasisY);
else
    newImage=[];
end






