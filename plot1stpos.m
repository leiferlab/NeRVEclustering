
function t_fh = plot1stpos(vtraj,traj,sptpara,options)

% INPUTS:
% vtraj = vector of traj numbers from spt output traj
%
% traj = structure of traj from spt
%
% sptpara = optional movie parameters
% sptpara.f_rate = frame rate of movie in hz
% sptpara.size = dimensions of movie in pixels ( [width height] )
% sptpara.pixl = pixel size in um/pixel
%
% options = optional structure of options
% options.overlay = if set to 1 plot is an overlay plot (no tick or axis
%                   labels, otherwise a normal graph is made
% options.lspec = set to a valid axis linespec in cell array form otherwise
%                 options.lspec.prop =
%                 {'Marker','MarkerSize','Line','Linewidth'}
%                 options.lspec.val = {'o',2,'-',1}
%                 is used.  
% options.color = set to 0 to plot trajs in color based on time, otherwise
%                 plots color based on traj id
%
% options.start = optional parameter frame number to start plotting (if not
%                 specfied start = 1
%
% options.stop = optional parameter frame number to stop plotting (if not 
%                specfied stop = inf
%
% OUTPUTS:
% t_fh = file handle of the figure created
%
% Dependencies: color_line
%
% Author: Colin Ingram
% Created: Aug 2009
% Version: 1.11
%
% Revisions:
% Version     Date  Author   Description
% 1.00   2009.08.26 CJI  Initial File
% 1.10   2010.04.01 SB   added option to color traj based on id
% 1.11   2010.04.05 CJI  made id color consistent with options format
%% Default Options

% default spt options
lim = 0;
frate = 1;
pixl = 1;

% default plot options
overlay = 0;
prop = {'Marker','MarkerSize','Line','Linewidth'};
val = {'.',6,'-',1};
color = 1;

% default times
st = 1;
sp = inf;


%% plots

% figure properties
fprop.NumberTitle = 'on';
fprop.PaperOrientation = 'landscape';
fprop.PaperPosition = [0 0.25 11 8];
fprop.PaperSize = [11 8.5];
fprop.Units = 'normalized';
fprop.Position = [.1 0 .9 .9];
t_fh = figure(fprop);
set(t_fh,'DefaultTextFontName','Arial');
set(t_fh,'DefaultAxesFontName','Arial');

% color properties (only used in identity plotting
colaxt = 'kgrbmc';
ntraj = size(vtraj,2);
ncolax = ceil(ntraj/size(colaxt,2));
colax = repmat(colaxt,1,ncolax);
    
for i = vtraj
    x = traj(i).pos(1,1);
    y = traj(i).pos(1,2);
    h = plot(x,y,'-','color',colax(i));
     set(h,prop,val);
    hold on
end
hold off    

% color bar
if color == 0
    colorbar
end


% axis properties
aprop.Box = 'on';
aprop.DataAspectRatio = [1 1 1];
aprop.YDir = 'reverse';

if lim == 1
    aprop.YLim = ylim;
    aprop.XLim = xlim;
end

if overlay == 1 
    aprop.YTickLabel = '';
    aprop.XTickLabel = '';
end

set(gca,aprop)
 print(gcf, '-dpdf', 'all traj.pdf')