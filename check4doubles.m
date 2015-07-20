function repTerms=check4doubles(x)
%check for repeated NONZERO elements of a vector
x=x(x~=0 & ~isnan(x));
xU=unique(x);
xhist=histc(x,xU);
repTerms=xU(xhist>1);



