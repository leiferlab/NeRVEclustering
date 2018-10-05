function B = internal_energy(alpha,beta, gamma, nPoints)

% Penta diagonal matrix, one row:
b(1)=beta;
b(2)=-(alpha + 4*beta);
b(3)=(2*alpha + 6 *beta);
b(4)=b(2);
b(5)=b(1);

% Make the penta matrix (for every contour point)
A=b(1)*circshift(eye(nPoints),2);
A=A+b(2)*circshift(eye(nPoints),1);
A=A+b(3)*circshift(eye(nPoints),0);
A=A+b(4)*circshift(eye(nPoints),-1);
A=A+b(5)*circshift(eye(nPoints),-2);
 A(end-1:end,1:2) = 0;
 A(1:2, end-1:end) = 0;
 
%  A(1,1:3) = [alpha/2 -alpha alpha/2];
%  A(2, 1:3) = -2*[alpha/2 -alpha alpha/2];
%  A(99,end-2:end) = -2*[alpha/2 -alpha alpha/2];
%  A(100,1:3) = [alpha/2 -alpha alpha/2];
 
 A(2,1:4) = [-alpha-2*beta, 2*alpha+5*beta, -alpha-4*beta, beta];
A(1,1:3) = [beta, -2*beta, beta];
  A(nPoints-1,end-3:end) = fliplr(A(2,1:4));
  A(nPoints, end-2:end) = [beta, -2*beta, beta];
  

% Calculate the inverse
B=inv(A + gamma.* eye(nPoints));
end