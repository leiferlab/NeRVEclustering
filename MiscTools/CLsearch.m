function output=CLsearch(midIm,CLX,CLY,show,x)

if nargin<4;
    show=0;
end
if nargin<5
    x=[0 0];
end


output=-nansum(interp2(midIm,CLX+x(1),CLY+x(2)));
output=output+x*x'/2500;
if show
imagesc(midIm)
hold on
plot(CLX+x(1),CLY+x(2),'x')
%axis equal
hold off
drawnow
end

