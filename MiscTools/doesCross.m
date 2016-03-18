function [bool, x, y]=doesCross(CLin)


m=diff(CLin);
m=m(:,2)./m(:,1);

dxmat=bsxfun(@minus,CLin(:,1),CLin(:,1)');
dymat=bsxfun(@minus,CLin(:,2),CLin(:,2)');

dv=bsxfun(@times, m, dxmat(1:end-1,:))-dymat(1:end-1,:);
dv2=diff(sign(dv),[],2);
dv2(isnan(dv2))=0;
dv3=abs(dv2)>1;
dv3=dv3.*dv3';
bool=any(dv3(:));
if nargout>1
    [x , y]=find(triu(dv3,1));
end