function[E,r]=plotMatch(x,y,yrange)
if nargin==2
    yrange=length(y);
end

yl=length(y);
xl=length(x);
mid=round(.5*(xl+yl));
counter=1;
r=-yrange:yrange;
E=zeros(1,length(r));
%% old method, a bit slower
% for i=r
%     if i>0
%     xtemp=x(i+1:end,:);
%     ytemp=y(1:yl-i,:);
%     elseif i==0
%         xtemp=x;
%         ytemp=y;
%     elseif i<0
%        ytemp=y(-i+1:end,:);
%     xtemp=x(1:yl+i,:);     
%         
%     end
%     E(counter)=-sum(1./sqrt(sum((xtemp-ytemp).^2,2)));
%     counter=counter+1;
%     
%  end
%%
X=convn(x,rot90(y,2),'full');
E=-X(:,2);
r=round((1:(yl+xl))-(yl+xl)/2);




% 
% plot(x(:,1),x(:,2))
% hold on
% plot(y(:,1),y(:,2))
% hold off
% axis equal
% 
% pause(1)


    
end
