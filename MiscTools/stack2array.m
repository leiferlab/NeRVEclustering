function [stack_matrix] = stack2array(stack_name)
%STACK2ARRAY takes a stacked .tif file and converts it to a 3D MATLAB
%matrix, especially for "WormSegmentHessianID.m".


%Convert a image file to an array of information structures.
info = imfinfo(stack_name);

%Finds the 3D dimensions in pixels of the stack, assuming each image in the
%stack has the same dimensions.
stack_z_size = length(info);
row = info(1).Height;
col = info(1).Width;

%Creates a 3D matrix with dimensions of the image stack. 
stack_matrix=zeros(row,col,stack_z_size);

%For each 2D level of the matrix, input the color information of the
%corresponding image in the stack as a double array.
for k = 1:stack_z_size
        stack_matrix(:,:,k) = double(imread(stack_name, 'tif', 'Index', k));
end


end