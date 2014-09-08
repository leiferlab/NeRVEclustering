im=uipickfiles();
%%
hyper_stack=  stackLoad(im{1});
imsize=size(hyper_stack);
%hyper_stack=image_resize(hyper_stack,imsize(1),imsize(2),imsize(3)*3.45);


%Convert the 3D matrix representation of the image stack to a 3D matrix
%representation of the binary image stack with parameters set in fields of
%"options".

%%
options.minObjSize=100;
options.maxObjSize=1000;
options.minSphericity=0.8;
options.watershedFilter=0;
options.filterSize=[40,40,10];
options.noise=.5;
options.maxSplit=0;
options.scaleFactor=[1,1,6];
thresh1=.05;




%%
[neurons,filteredWorm] = WormSegmentHessian3d_rescale(hyper_stack, options);


%%
neuronsLabel=bwlabeln(neurons);