%   Options.SigmaFluid : The sigma of the gaussian smoothing kernel of the pixel
%                   velocity field / update field, this is a form of fluid
%                   regularization, (default 4)
%   Options.SigmaDiff : The sigma for smoothing the transformation field
%                   is not part of the orignal demon registration, this is 
%                   a form of diffusion regularization, (default 1)
%   Options.Interpolation : Linear (default) or Cubic.
%   Options.Alpha : Constant which reduces the influence of edges (and noise)
%                   and limits the update speed (default 4). 
%   Options.Similarity : Choose 'p' for single modality and 'm' for
%                   images of different modalities. (default autodetect)
%   Options.Registration: Rigid, Affine, NonRigid  
%   Options.MaxRef : Maximum number of grid refinements steps.
%   Options.Verbose: Display Debug information 0,1 or 2 
%    
%   % Read two greyscale images of Lena
vidFile=uipickfiles();
vidFile=vidFile{1};
VidObj= VideoReader(vidFile);
%%
Istatic=read(VidObj,1);
Imoving=read(VidObj,100);
Imoving=im2double(Imoving(400:700,400:700,1));

Istatic=im2double(Istatic(400:700,400:700,1));
%%


%   Imoving=imread('images/lenag1.png'); 
%   Istatic=imread('images/lenag3.png');
%   
  Options.Registration='NonRigid';
  Options.Similarity='p';
  Options.MaxRef=50;
Options.Verbose=0;
  % Register the images
  tic
  [Ireg,Bx,By,Fx,Fy] = register_images(Imoving,Istatic,Options);
toc
  % Show the registration result
  figure,
  subplot(2,2,1), imshow(Imoving); title('moving image');
  subplot(2,2,2), imshow(Istatic); title('static image');
  subplot(2,2,3), imshow(Ireg); title('registerd moving image');
  % Show also the static image transformed to the moving image
  Ireg2=movepixels(Istatic,Fx,Fy);
  subplot(2,2,4), imshow(Ireg2); title('registerd static image');

 % Show the transformation fields
  figure,
  subplot(2,2,1), imshow(Bx,[]); title('Backward Transf. in x direction');
  subplot(2,2,2), imshow(Fx,[]); title('Forward Transf. in x direction');
  subplot(2,2,3), imshow(By,[]); title('Backward Transf. in y direction');
  subplot(2,2,4), imshow(Fy,[]); title('Forward Transf. in y direction');

% Calculate strain tensors
  E = strain(Fx,Fy);
% Show the strain tensors
  figure,
  subplot(2,2,1), imshow(E(:,:,1,1),[]); title('Strain Tensors Exx');
  subplot(2,2,2), imshow(E(:,:,1,2),[]); title('Strain Tensors Exy');
  subplot(2,2,3), imshow(E(:,:,2,1),[]); title('Strain Tensors Eyx');
  subplot(2,2,4), imshow(E(:,:,2,2),[]); title('Strain Tensors Eyy');

