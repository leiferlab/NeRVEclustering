function[E,r]=plotMatch(x,y,yrange)
if nargin==2
    yrange=length(y);
end

yl=length(y);
xl=length(x);
r=-yrange:yrange;
% w=gausswin(length(x),3);
% w=w/sum(w);
x=bsxfun(@minus,x,mean(y));
y=bsxfun(@minus,y,mean(y));

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
xref=convn(ones(size(x)),ones(size(x)),'full');
 X=X./xref;
E=-X(:,2);
r=round((1:(yl+xl))-(yl+xl)/2);





    
end
