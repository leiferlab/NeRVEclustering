function Fline = external_energy(I, sigma)

%Make derivative kernals
[x,y]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));

DGausx =-(x./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
DGausy = -(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));

% Calculate the forces from the image

 Fline(:,:,1) = -imfilter(I, DGausx, 'conv', 'symmetric')*sigma;
 Fline(:,:,2) = -imfilter(I, DGausy, 'conv', 'symmetric')*sigma;
 
end