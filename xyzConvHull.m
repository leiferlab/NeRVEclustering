function [ Iout ] = xyzConvHull( I ,dims)
%UNTITLED4 use BWfill in x y and z, or selected dimensions.
I=imdilate(I,true(3,3,3));
if nargin==1
    dims=1:3;
end

for runs=1:2
    if any(dims==1)
        for imSlice=1:size(I,1);
            I(imSlice,:,:)=imfill(squeeze(I(imSlice,:,:)),'holes');
        end
    end
    if any(dims==3)
        for imSlice=1:size(I,3);
            I(:,:,imSlice)=imfill(I(:,:,imSlice),'holes');
        end
    end
    if any(dims==2)
        for imSlice=1:size(I,2);
            I(:,imSlice,:)=imfill(squeeze(I(:,imSlice,:)),'holes');
        end
    end
    
end



Iout=imerode(I,true(3,3,3));