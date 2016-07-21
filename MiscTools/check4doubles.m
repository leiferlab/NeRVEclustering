function repTerms=check4doubles(x)
%returns non-zero and non-nan elements of an array that appear more than
%once
%remove zeros and nans
x=x(x~=0 & ~isnan(x));
%find unique terms and their indices
[x_unique,~,ib]=unique(x);
% count how many of each index appears
xhist=accumarray(ib(:),ones(length(ib),1));
%get those values
repTerms=x_unique(xhist>1);



