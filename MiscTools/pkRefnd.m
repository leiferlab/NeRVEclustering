% refines the pixel level accurancy peaks from pkfnd by applying a
% threshold to the integrated intensity of a region surrounding the feature
% centers located within half the region length of the edges will not be
% consider.
function out=pkRefnd(frame,pkpeaks,boxr,thresh)
% INPUTS:
% frame: image to process, particle should be bright spots on dark 
%        background with little noise ofren an bandpass filtered 
%        brightfield image (fbps.m, fflt.m or bpass.m) or a nice 
%        fluorescent image
% pkpeaks: list of particle centers to pixel level accuracy (pkfnd.m
%          output [row,column])
% boxr: length of box defineing particle region
% th:  the minimum brightness of region containing particle

% OUTPUTS:  a N x 3 array containing, [row,column,intensity] 
%         coordinates of local maxima
%           out(:,1) are the x-coordinates of the maxima
%           out(:,2) are the y-coordinates of the maxima
%           out(:,3) are the integrated intensity of the region containing
%           the maxima
%
% Dependencies: None
%
% Author: Colin Ingram
% Created: May 2010
% Version: 0.1
%
% Revisions:
% Version  Date       Author   Description
% 0.01     2010.04.11 CJI      Initial Creation 
% 0.02     2010.05.04 CJI      fixed edge protection

%% Prep and Inputs
sz = floor(boxr/2);
[nr,nc]=size(frame);

% remove peaks to close to the edges
ind = pkpeaks;
ind=ind(ind(:,2) > sz & ind(:,2) < nr-sz,:);
ind=ind(ind(:,1) > sz & ind(:,1) < nc-sz,:);
[nind,~] = size(ind);

%% process
for j = 1:nind
    x = ind(j,1);
    y = ind(j,2);
    ind(j,3) = sum(sum(frame(y-sz:y+sz,x-sz:x+sz)));
end

out = ind(ind(:,3)>thresh,:);
