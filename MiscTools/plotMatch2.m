function[E,r]=plotMatch2(x,y,yrange)
%finds offset between two 2D curves so that points align better, outputs an
%energy E for each offset r
if nargin==2
    yrange=length(y);
end

yl=length(y);
xl=length(x);
r=-yrange:yrange;
% w=gausswin(length(x),3);
% w=w/sum(w);
dmat=pdist2(y,x);
[xmat,ymat]=meshgrid(1:xl,1:yl);
diffmat=xmat-ymat;
E=accumarray(diffmat(:)+xl,dmat(:),[],@mean);
r=-(xl-1):(xl-1);





    
end
