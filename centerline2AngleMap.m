function [dangle, f]=centerline2AngleMap(centerline)

%function takes centerline data from ashley an N by 2 by Time matrix and
%produces an N by Time matrix of angles
%%
if ndims(centerline)==3
tangentVector=diff(centerline,[],1);
tangentVector=bsxfun(@rdivide,tangentVector,...
    hypot(tangentVector(:,1,:),tangentVector(:,2,:)));
%deltaT=diff((tangentVector),[],1);
deltaT=tangentVector;
angles=atan2(deltaT(:,1,:),deltaT(:,2,:));
angles=unwrap(angles);
dangle=diff(angles,[],1);
dangle=squeeze(dangle);
elseif ismatrix(centerline)
    dangle=centerline;
end
x=1:size(dangle,1);

a=zeros(size(x));b=a;c=a;
progressbar(0)
for i=1:size(dangle,2);
try    
    progressbar(i/size(dangle,2));
    f=fit(x',dangle(:,i),'sin1');
    a(i)=f.a1;b(i)=f.b1;c(i)=f.c1;
catch
end

end
c=unwrap(c);
v=diff(smooth(c,200));
f=struct('a',a,'b',b,'c',c,'v',v);