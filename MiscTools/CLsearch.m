function output=CLsearch(midIm,CLX,CLY,show,x)
%interpolates the values of midIm at points [CLX,CLY], given an offset of
%x. Inverts the sign of the output
if nargin<4;
    show=0;
end
if nargin<5
    x=[0 0];
end


output=-nansum(interp2(midIm,CLX+x(1),CLY+x(2)));
output=output+x*x'/2500; %small energetic penalty for large offsets
if show
imagesc(midIm)
hold on
plot(CLX+x(1),CLY+x(2),'blackx')
%axis equal
hold off
drawnow
end

