function output=CLsearch_2(midIm,CLX,CLY,show,x)
%%same as CLsearch, which interpolates the values of in image at certain
%%coordinates and returns the sum of the interpolated values. This version
%%speeds things up by not using interp funciton, instead just using
%%indexing

if nargin<4;
    show=0;
end
if nargin<5
    x=[0 0];
end
s_ind=[round(CLX) round(CLY)];

[row,col]=size(midIm);


        in_image_flag = all((s_ind>0)&(bsxfun(@le, s_ind, [col, row])), 2);
        gradient_ind = sub2ind_nocheck([row, col], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1));
output=-nansum(midIm(gradient_ind));
output=output+x*x'/2500;
if show
imagesc(midIm)
hold on
plot(CLX+x(1),CLY+x(2),'blackx')
%axis equal
hold off
drawnow
end

