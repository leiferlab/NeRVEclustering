function repTerms=check4doubles(x)
%returns non-zero and non-nan elements of an array that appear more than
%once
x=x(x~=0 & ~isnan(x));
xU=unique(x);
xhist=histc(x,xU);
repTerms=xU(xhist>1);



