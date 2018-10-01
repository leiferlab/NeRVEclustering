function mask=FourierMaskHexagon(hexlength,img,holesize)
% Generate fourier mask for a hexagon lattice.(already fftshift)
% The length of size of triangle is length.

% Input:
%   length: the length of triangle side(pixel).
%   img: the image contains hexagon lattice. To get image size.
%   holesize: the size of hole in the mask
% Output:
%   mask: mask show where the fourier frequency peak is.( use fftshift to
%       center the maximum value)

    if nargin<3
        holesize=7;
    end
    
    % Get the size of image
    [image_width,image_height]=size(img);
    % start at the center of image,x is the first dimension(height) y is the
    % second dimension(width)
    start_y=ceil( (1+image_width)/2 );
    start_x=ceil( (1+image_height)/2 );
    
    % Get the period for Fourier frequency
    % Got by observed from standard lattice. Make sense if considering the
    % lattice period.
    % frequency considering x and y direction
    Fourier_y=image_width/hexlength;
    Fourier_x=image_height/(hexlength*sqrt(3));
    
    
    mask=zeros(image_width,image_height);

    
  % Need to do it in four direction
  
  %First direction: up right
    % row index is the order of row in the lattice.
    % row is the real coordinate of the row
    rowindex=1;
    row=start_x;

    while row<=image_height
        % odd row and even row is not exact same.
        if mod(rowindex,2)
            column=start_y;
        else
         % add the minimum frequency of x direction.
            column=start_y+Fourier_y;
        end
    
        while column<=image_width
            mask(floor(row),floor(column))=1;
            % times 2 because the properties of hexagonal lattice. (not
            % exact orthogonal but has two pattern repeat.
            column=column+2*Fourier_y;
        end
        % next row, add sqrt(3)/2 * length of triangle
        row=row+Fourier_x;   
        rowindex=rowindex+1;
    end
    
    %Second direction: up left
    % row index is the order of row in the lattice.
    % row is the real coordinate of the row
    rowindex=1;
    row=start_x;

    while row<=image_height
        % odd row and even row is not exact same.
        if mod(rowindex,2)
            column=start_y;
        else
         % add the minimum frequency of x direction.
            column=start_y-Fourier_y;
        end
    
        while column>=1
            mask(floor(row),floor(column))=1;
            % times 2 because the properties of hexagonal lattice. (not
            % exact orthogonal but has two pattern repeat.
            column=column-2*Fourier_y;
        end
        % next row, add sqrt(3)/2 * length of triangle
        row=row+Fourier_x;   
        rowindex=rowindex+1;
    end
    
  %Third direction: down right
    % row index is the order of row in the lattice.
    % row is the real coordinate of the row
    rowindex=1;
    row=start_x;

    while row>=1
        % odd row and even row is not exact same.
        if mod(rowindex,2)
            column=start_y;
        else
         % add the minimum frequency of x direction.
            column=start_y+Fourier_y;
        end
    
        while column<=image_width
            mask(floor(row),floor(column))=1;
            % times 2 because the properties of hexagonal lattice. (not
            % exact orthogonal but has two pattern repeat.
            column=column+2*Fourier_y;
        end
        % next row, add sqrt(3)/2 * length of triangle
        row=row-Fourier_x;   
        rowindex=rowindex+1;
    end  
    
  %Fourth direction: down left
    % row index is the order of row in the lattice.
    % row is the real coordinate of the row
    rowindex=1;
    row=start_x;

    while row>=1
        % odd row and even row is not exact same.
        if mod(rowindex,2)
            column=start_y;
        else
         % add the minimum frequency of x direction.
            column=start_y+Fourier_y;
        end
    
        while column>=1
            mask(floor(row),floor(column))=1;
            % times 2 because the properties of hexagonal lattice. (not
            % exact orthogonal but has two pattern repeat.
            column=column-2*Fourier_y;
        end
        % next row, add sqrt(3)/2 * length of triangle
        row=row-Fourier_x;   
        rowindex=rowindex+1;
    end    
   
    % dilate surrounding area.
    se=strel('square',holesize);
    
    mask=imdilate(mask,se);
end


