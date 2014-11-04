

p = 'F:\20141029\chip_sp_gc6rfp_1';

files = dir([p filesep '*.tif']);

image = imread([p filesep files(1).name], 'tif');


imagesc(image);
rect1 = getrect;
% rect2 = getrect;
% rect3 = getrect;

rect1 = round(rect1);
% rect2 = round(rect2);
% rect3 = round(rect3);


ROIval1 = zeros(1, length(files));
% ROIval2 = zeros(1, length(files));
% ROIval3 = zeros(1, length(files));

for iImage = 1:length(files)
    
   temp =  imread([p filesep files(iImage).name], 'tif');
   temp1 = temp(rect1(2):rect1(2)+rect1(4), rect1(1):rect1(1)+rect1(3));
   
   ROIval1(iImage) = mean(temp1(:));
   
%    
%   
%    temp2 = temp(rect2(2):rect2(2)+rect2(4), rect2(1):rect2(1)+rect2(3));
%    
%    ROIval2(iImage) = mean(temp2(:));
%    
%    temp3 = temp(rect3(2):rect3(2)+rect3(4), rect3(1):rect3(1)+rect3(3));
%    
%    ROIval3(iImage) = mean(temp3(:));
    
end