function [circleEnergy,projections]=circlePlane(X,theta,phi)
%function is used to find best plane to fit circular data. to be used with
%an fminserach


n=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

T1=normr(cross(n, -circshift(n, [0 1])));
T2=normr(cross(n,T1));
projectionX=X*T1';
projectionY=X*T2';
projectionX=projectionX-mean(projectionX);
projectionY=projectionY-mean(projectionY);

R=sqrt(projectionX.^2+projectionY.^2);
%  scatter(projectionX, projectionY)
%  drawnow
circleEnergy=std(R);
projections=[projectionX, projectionY];