im=uipickfiles();
%%
hyper_stack=  stackLoad(im{1});
imsize=size(hyper_stack);
hyper_stack=image_resize(hyper_stack,imsize(1),imsize(2),imsize(3)*3.45);


%Convert the 3D matrix representation of the image stack to a 3D matrix
%representation of the binary image stack with parameters set in fields of
%"options".

%%
options.minObjSize=100;
options.maxObjSize=200;
options.minSphericity=0.9;
options.watershedFilter=1;
options.filterSize=[30,30,15];
options.noise=.5;

%%
neurons = WormSegmentHessian3D(hyper_stack, options);