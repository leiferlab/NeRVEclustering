function Pout=coordinateTransform3d(P,X,Y,Z)
%%%% 
% coordinateTransform3d takes a Nx3 point and 3D lookup tables for X, Y and
% Z inorder to find the new Nx3 coordinates. 
imSize=size(X);
Plin=round(P);
Plin(Plin<1)=1;
Plin(Plin(:,1)>=imSize(1),1)=imSize(1);
Plin(Plin(:,2)>=imSize(2),2)=imSize(2);
Plin(Plin(:,3)>=imSize(3),3)=imSize(3);
Plin=sub2ind(imSize,Plin(:,1),Plin(:,2),Plin(:,3));


Px=X(Plin);
Py=Y(Plin);
Pz=Z(Plin);

% Px=interp3(X,P(:,2),P(:,1),P(:,3),'nearest*');
% Py=interp3(Y,P(:,2),P(:,1),P(:,3),'nearest*');
% Pz=interp3(Z,P(:,2),P(:,1),P(:,3),'nearest*');
Pout=[Px Py Pz];
