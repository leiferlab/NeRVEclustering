function [Iout,lookup]=imcrop3d(I,rect)
%quick crop of 3d volume from bounding rect from regionprops
rect=round(rect);
rect(4:end)=rect(4:end)-1;
Iout=I(rect(2)+(0:rect(5)),rect(1)+(0:rect(4)),rect(3)+(0:rect(6)));

if nargout==2
    [x, y, z]=ndgrid(rect(1)+(0:rect(4)),rect(2)+(0:rect(5)),rect(3)+(0:rect(6)));
lookup=sub2ind(size(I),y(:),x(:),z(:));
end
