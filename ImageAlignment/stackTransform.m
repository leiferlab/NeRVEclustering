function Vout=stackTransform(V,tformcell)

% stackTransform takes a MxNxL image V and a 1XL cell array of tforms and
% performs the tform on each of the image slices of V. The output image has
% the same size as the input image.
imsize=size(V);
Vout=V;
V(isnan(V))=0;
 R = imref2d(size(V(:,:,1))) ;
for iSlice=1:imsize(3)
    Vout(:,:,iSlice)=imwarp(V(:,:,iSlice),R,tformcell{iSlice},...
    'nearest','OutputView',R);
    
end


