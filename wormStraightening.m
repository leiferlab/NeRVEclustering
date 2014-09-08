function [newImage,newBasisX,newBasisY]=wormStraightening(CL,imageIn,windowSize);

windowSize=20;
windowSpan=-windowSize:1:windowSize;
[~,perpSlopes]=gradient(CL,5);
perpSlopes=normr([-perpSlopes(:,2),perpSlopes(:,1)]);


newBasisX=perpSlopes(:,1)*windowSpan;
newBasisY=perpSlopes(:,2)*windowSpan;

newBasisX=bsxfun(@plus,newBasisX,CL(:,1));
newBasisY=bsxfun(@plus,newBasisY,CL(:,2));

% imagesc(imageIn)
% hold on
% plot(CL(:,2),CL(:,1),'b');
% plot(newBasisY',newBasisX','black')


newImage=interp2(sum(imageIn,3)',newBasisX,newBasisY);






