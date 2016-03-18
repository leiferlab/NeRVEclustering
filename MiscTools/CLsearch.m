function output=CLsearch(midIm,CLX,CLY,show)

if nargin==3;
    show=0;
end

output=-nansum(interp2(midIm,CLX,CLY));
if show
imagesc(midIm)
hold on
plot(CLX,CLY,'x')
%axis equal
hold off
drawnow
end

