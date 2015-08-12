function varargout=  stackLoad(stack_name, stack_z_size,Z_scale)
% takes a tiff stack and loads them into 3D matrices. Will take the entire
% stack and load them in chuncks of stack_z_size, if the stack is larger.
% Use multiple outputs if you expect multiple stacks because of this.
% default stack_z_size is size of entire stack. If there are no outputs,
% will display matrix with sliceBrowser.
[row, col] = size(imread(stack_name, 'tif', 1));
if nargin==1;

totalStack=length(imfinfo(stack_name));
stack_z_size =totalStack;
Z_scale=1;         %.6   %scaling between distance in the slide to distance moved in the image
end

if nargin==2
Z_scale=1;         %.6   %scaling between distance in the slide to distance moved in the image
if isempty(stack_z_size)
    totalStack=length(imfinfo(stack_name));
    stack_z_size=totalStack;
else
    totalStack=stack_z_size;
end
end

for h=1:min(floor(totalStack/stack_z_size), nargout)
    hyper_stack=zeros(row,col,stack_z_size);
    for k = 1:stack_z_size
        hyper_stack(:,:,k) = double(imread(stack_name, 'tif',...
            (h-1)*stack_z_size+k));
    end
    if Z_scale~=1
    hyper_stack=image_resize(hyper_stack,size(hyper_stack,1),...
        size(hyper_stack,2),round(Z_scale*size(hyper_stack,3))); 
    end
    
    if nargout==0
    SliceBrowser(hyper_stack)
    end
    varargout{h}=hyper_stack;
    
end
