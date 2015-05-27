function MI=mutualInformation(varargin)
if nargin>=2
    x=varargin{1};y=varargin{2};
    if nargin>=3
        e1=varargin{3};
        if nargin>=4
            e2=varargin{4};
        else e2=10;
        end
    else e1=10;
    end
    
    
    
for i=1:size(y,1);
    ytemp=y(i,:);
    xtemp=x;
%     bad=isnan(ytemp(:)+xtemp(:));
%     ytemp=detrend(ytemp(~bad));
%     xtemp=detrend(xtemp(~bad));
dmat=histcn([xtemp(:),ytemp(:)],e1,e2);
dmat=dmat/nansum(dmat(:));
MI(i)= nansum(nansum(dmat.*log2(dmat./(bsxfun(@times,nansum(dmat,1),nansum(dmat,2))))));

end
elseif nargin==1
    dmat=varargin{1};
    MI= nansum(nansum(dmat.*log2(dmat./(bsxfun(@times,nansum(dmat,1),nansum(dmat,2))))));

end
