function tiffwrite(filename,imageData,softwareTag,showProgressBar);
% Inputs:
%     filename, a string with relative or full path of where to write
%         the tif image
%     image, a m x n x k array of values to be written as a tiff
%     softwareTag (optional), a string to be used for the software tif tag
%     showProgressBar (optional), a boolean value for showing a progress
%             bar, default is false and a value of true requires progressbar
%             from matlab central
%
% Dependencies: progressbar from mathworks file exchange;
% Author: Benjamin Bratton, UW Madison, bpbratton@gmail.com
% Version: 0.1 20110330
% Comments: BPB Compared writing images with Tiff and with imwrite and Tiff
%   was significantly faster for large images, especially with many many
%   planes
%

t = Tiff(filename,'w');
tagstruct.ImageLength = size(imageData,1);
tagstruct.ImageWidth = size(imageData,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;

switch class(imageData);
    case 'double'
        warning('tiffwrite:imgClass:double',...
            'unsupported image input class of ''double'' was recast as single');
        imageData = single(imageData);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'single'
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'int8'
        tagstruct.BitsPerSample = 8;
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
    case 'int16'
        tagstruct.BitsPerSample = 16;
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
    case 'uint8'
        tagstruct.BitsPerSample = 8;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case 'uint16'
        tagstruct.BitsPerSample = 16;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case 'int32'
        warning('tiffwrite:imgClass:int32',...
            'unsupported image input class of ''int32'' was recast as single');
        imageData = single(imageData);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'uint32'
        warning('tiffwrite:imgClass:int32',...
            'unsupported image input class of ''uint32'' was recast as single');
        imageData = single(imageData);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'int64'
        warning('tiffwrite:imgClass:int32',...
            'unsupported image input class of ''int64'' was recast as single');
        imageData = single(imageData);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'uint64'
        warning('tiffwrite:imgClass:int32',...
            'unsupported image input class of ''uint64'' was recast as single');
        imageData = single(imageData);
        tagstruct.BitsPerSample = 32;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'logical'
        tagstruct.BitsPerSample = 1;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
end

if nargin>2
    if not(isempty(softwareTag))
        tagstruct.Software = softwareTag;
    end
end

if nargin<4
    showProgressBar = false;
end

if showProgressBar
    progressbar('Writing tif');
end

t.setTag(tagstruct);
t.write(imageData(:,:,1));
nPlanes = size(imageData,3);

for ind  = 2:nPlanes;
    if showProgressBar
        progressbar((ind-1)/nPlanes);
    end
    t.writeDirectory();
    t.setTag(tagstruct);
    t.write(imageData(:,:,ind));
end
t.close();

if showProgressBar
    progressbar(1);
end

