masterFiducials=rand(65,3)-.5;
masterFiducials=bsxfun(@minus,masterFiducials,mean(masterFiducials));
currentFiducials=masterFiducials+(rand(size(masterFiducials))-.5)*0.2+1;
%PointsMoving=fliplr(PointsMoving);
currentFiducials=circshift(currentFiducials,[0 1]);
%%

DT= delaunayTriangulation(masterFiducials);

CMatrix=false(size(masterFiducials,1));
for i=1:length(masterFiducials)
    [x,y]=find(DT.ConnectivityList==i);
    rows=DT.ConnectivityList(x,:);
    rows=unique(rows(:));
    idx=sub2ind(size(CMatrix),repmat(i,1,length(rows))',rows);
    CMatrix(idx)=true;
    
end
CMatrix=CMatrix.*~eye(length(CMatrix));
%%
DmatX0=bsxfun(@minus, masterFiducials(:,1),masterFiducials(:,1)');
DmatY0=bsxfun(@minus, masterFiducials(:,2),masterFiducials(:,2)');
DmatZ0=bsxfun(@minus, masterFiducials(:,3),masterFiducials(:,3)');

Dmat0=sqrt(DmatX0.^2+DmatY0.^2+DmatZ0.^2);
CMatrix(Dmat<100)=0;

%%
PointsMovingNew=currentFiducials;
Lmatrix=-CMatrix;
Lmatrix(eye(length(CMatrix))>0)=-sum(Lmatrix);
  %%

      tform = makeAffine3d(masterFiducials,PointsMovingNew );

           [masterFiducialsMoving(:,1),masterFiducialsMoving(:,2),masterFiducialsMoving(:,3)]...
               =transformPointsForward(tform,masterFiducials(:,1),...
               masterFiducials(:,2),masterFiducials(:,3));   
           
           
           
     

%%
 close all
 gamma=1;
 kappa=.001;
 controlIdx=1:2:length(masterFiducialsMoving);
 noise=.00;
 maxIt=5000;
for iteration=1:maxIt;
DmatX=bsxfun(@minus, masterFiducialsMoving(:,1),masterFiducialsMoving(:,1)');
DmatY=bsxfun(@minus, masterFiducialsMoving(:,2),masterFiducialsMoving(:,2)');
DmatZ=bsxfun(@minus, masterFiducialsMoving(:,3),masterFiducialsMoving(:,3)');
Dmat=sqrt(DmatX.^2+DmatY.^2+DmatZ.^2);
% Hcell=[];
% for i=1:length(DmatX);
%     for j=1:length(DmatX)
%     Hij=[DmatX(i,j) DmatY(i,j) DmatZ(i,j)]'* [DmatX(i,j) DmatY(i,j) DmatZ(i,j)];
%     Hcell{i,j}=-Hij/Dmat(i,j).^2;
%     endxx
% end
% 

% H=cell2mat(Hcell);
% H(isnan(H))=0;
% H=H*.01;
% 
%     PointsMovingNewTemp=PointsMovingNew';
%     
%     PointsMovingNewTemp=PointsMovingNewTemp(:);
%     
%     f=(speye(size(H))-H);
% [FL,FU]=lu(f);
%     
%     
%      PointsMovingNewTemp = FU\(FL \...
%         (PointsMovingNewTemp));
%          F=reshape(PointsMovingNewTemp,3,[])'-PointsMovingNew;
% 
%     PointsMovingNew=reshape(PointsMovingNewTemp,3,[])';
   % PointsMovingNew=Points+new_deltaPosition;
%     Fmatrix=Lmatrix.*eye(length(Lmatrix))+Lmatrix;
%     H

  %  PointsMovingNew(controlIdx,:)=Points(controlIdx,:);
  %  V=sum(diag((Lmatrix*deltaPosition)'*deltaPosition));
    V=sum(sum(triu(Dmat-Dmat0).^2))
    S=1-Dmat0./Dmat;
    
    S(isnan(S))=0;
%     H=-kappa*(S*Lmatrix);
%     
%     F=Lmatrix*PointsMovingNew./Dmat0.*(Dmat-Dmat0).*Cmatrix;
% 
%         f=(speye(size(H))*gamma+H);
% [FL,FU]=lu(f);
%     
%     
%      PointsMovingNewTemp = ...%FU\(FL \...
%         (f\(gamma*PointsMovingNew));
%          F=(PointsMovingNewTemp-PointsMovingNew);
%          
% 
%     F=bsxfun(@minus,F,mean(F));
%    PointsMovingNew(1:10,:)= PointsMovingNew(1:10,:)+F(1:10,:);
% 

   Fx=sum(S.*(DmatX).*CMatrix);
   Fy=sum(S.*(DmatY).*CMatrix);
   Fz=sum(S.*(DmatZ).*CMatrix);
   
% deltaX=DmatX-DmatX0;
% deltaY=DmatY-DmatY0;
% % deltaZ=DmatZ-DmatZ0;
% % 
% deltaS=(DmatX.*deltaX + DmatY.*deltaY+DmatZ.*deltaZ);
% deltaS=.5*sign(deltaS)./sqrt(abs(Dmat));
% deltaS=1-Dmat0./Dmat;
% Fx=deltaS.* (DmatX)./Dmat.*CMatrix;
% Fy=deltaS.* ( DmatY)./Dmat.*CMatrix;
% Fz=deltaS.* (DmatZ)./Dmat.*CMatrix;
% 
% Fx=-nansum(Fx)/2;
% Fy=-nansum(Fy)/2;
% Fz=-nansum(Fz)/2;
% % 

% % 

F=[Fx' Fy' Fz']*kappa;
Ffiducials=zeros(size(F));
F=F+(rand(size(F))-.5)*noise*(1-iteration/maxIt);
  masterFiducialsMoving=masterFiducialsMoving+F;
  
  
 masterFiducialsMoving=bsxfun(@minus,masterFiducialsMoving,mean(masterFiducialsMoving));
  masterFiducialsMoving=bsxfun(@plus,masterFiducialsMoving,mean(PointsMovingNew));
Ffiducials(controlIdx,:)=PointsMovingNew(controlIdx,:)-masterFiducialsMoving(controlIdx,:);
masterFiducialsMoving=masterFiducialsMoving+Ffiducials*kappa*10;
  
  %masterFiducialsMoving(controlIdx,:)=PointsMovingNew(controlIdx,:);
  if mod(iteration,500)==0
delete(gca)
scatter3(PointsMovingNew(:,1),PointsMovingNew(:,2),PointsMovingNew(:,3));
hold on
scatter3(masterFiducialsMoving(:,1),masterFiducialsMoving(:,2),masterFiducialsMoving(:,3),'rx');
quiver3(masterFiducialsMoving(:,1),masterFiducialsMoving(:,2),masterFiducialsMoving(:,3),...
    F(:,1),F(:,2),F(:,2));
%axis([0 1 0 1 0 1]-.5);
drawnow
end
end



