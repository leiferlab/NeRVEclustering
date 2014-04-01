% spt_fndpos will calculate the trajectories of single particles in a stack
% tiff movie.  It is made to be used with spt.m but can be used on its on
% as well. Note: Some error checking is done in spt.m wrapper so only use
% this if you know what you are doing :-)

function [im, poslist] = spt_fndpos(sptpara,im,warn)

% INPUTS:
% sptpara.thresh = particle threshhold: use by pkfind.m to find the coarse
%                  pixle locations of particles
%
% sptpara.IntTh = [OPTIONAL] enables the use of refined threshholding.
%                 coarse locations will only be consider if the integrated
%                 intensity of the region within sptpara.boxr is greater
%                 than sptpara.IntTh. [default sptpara.IntTh=0]
%
% sptpara.dia = typical particle diameter: this is the diameter of the 
%               of a typical particle in pixels (MUST BE ODD).  
%               This parameter is used in several of the fuctions called 
%               by spt_calc. See those files for additional information. 
%
%               bpass.m: lobject = sptpara.dia+2 defines characteristics 
%                        length scale of objects to keep when smoothing
%               pkfind.m: sz = sptpara.dia+2 defines a radius (sz/2). If
%                         multiple peaks are found with in this radius only
%                         one is kept
%               centrd: sz = sptpara.dia+4 defines the length of window
%                            which should contain the entire particle but 
%                            no others. The centroid calculation is done in
%                            this window.  Also a particle within this
%                            distance from the edge is eliminated.
%
% sptpara.max_disp = [OPTIONAL] A warning will be
%                    issued if 2*max_disp > average particle distance.
%                    Since the data begains to become questionable.  As
%                    max_disp approaches the average particle distance the
%                    amount of time needed to calcuate traj ballons and
%                    sometimes track will fail. 
%
% sptpara.file = stacked tiff filename: the stacked tiff filename (string)
% sptpara.boxr = size of fit box for centfind, must be odd.
% sptpara.fitmethod = set to zero to use gaussian fit method, set to one to
%                     use centroid method
% im = structured array containing the raw image as
%   im.rawImg(height,width,nPlanes)
%
% warn[optional] = set to 0 to turn off warnings, useful for
%                  spt_interactive
%
% OUTPUTS:
% im = input im array plus im.planeAttr(:).* with entries for each plane     
%      particlepos: an array [xi,yi] of particle positions
%      nparticle: the number of particles found
%      iparticledis: a vector containing the mean interparticle distance(1)
%                    and its standard deviations(2)
% im also contains im.filterImg(height,width,nPlanes), the band pass
%       filtered image
%      
% trajlist = an array contain list of traj found [xi, yi, idf, idt], xi
%            is is the xi, yi is the coordinate of the particle which is 
%            located in frame idf and belongs to trajectory idt.
%
% Dependencies: bpass.m
%               pkfind.m
%               pkRefnd.m
%               centfind.m
%
% Author: Colin Ingram
% Created: June 2007
% Version: 1.01
%
% Revisions:
% Version     Date  Author   Description
% 0.01    2007.06    CJI: Initial Creation 
% 0.02    2008.10.28 CJI: moved filename into sptpara structure
% 1.00    2008.11.24 CJI: minor bug fixes, added interdistance calculation
%                         and warnings, and im outputs
% 1.01    2009.09.30 RMD: fixed diameter bug
% 1.50    2010.02.28 SB:  handles frames without any particles correctly
% 1.60    2010.03.12 CJI: integrated intensity threshold added
% 2.00    2010.03.13 CJI: split file into spt_track and spt_fndpos
% 2.01    2010.05.04 CJI: added bpsz stuff to centfind
% 2.02    2010.07.16 BPB: changed structure of im
% 2.02    2011.03.21 BPB: updated documentation about new im structure

%%

% Handle Inputs
thresh = sptpara.thresh;
dia = sptpara.dia;
bpsz = dia+2;
filename = sptpara.file;
opt.method = sptpara.fitmethod;
boxr = sptpara.boxr;

if ~exist('warn','var')
    warn = 1;
end

if isfield(sptpara,'IntTh')
    intth = sptpara.IntTh;
else
    intth = 0;
end

if isfield(sptpara,'max_disp')
    max_disp = sptpara.max_disp;
end

% im size
imcnt = size(im.rawImg,3);

% preallocate matrix for filtered images
im.filterImg = zeros(size(im.rawImg),'single');

% Find Particle Positions
t_pos = [];

% Suppress warnings from pkfind
wstate1 = warning('query','pkfnd:nothingAboveThresh');
warning('off','pkfnd:nothingAboveThresh');
for i = 1:imcnt
    
    % Clear outputs which maybe set (spt_interactive)
    im.planeAttr(i).iparticledis = [];
    im.planeAttr(i).particlepos = [];
    im.planeAttr(i).nparticle = [];
    im.planeAttr(i).particleintthr = [];
    im.planeAttr(i).particlepxl = [];
    
    
    % filter image and find coarse centers
    % (CJI+RMD 2009.09.30) bpass and pkfind request diameter slightly larger than
    % object diameter
    im.filterImg(:,:,i) = bpass(im.rawImg(:,:,i),1,dia+2,0,'single');
    t_pkfind = pkfnd(im.filterImg(:,:,i),thresh,bpsz);
    
    % Warn on empty frames
    if isempty(t_pkfind)
        
        % set nparticle to zero
        im.planeAttr(i).nparticle = 0;
        
        % warn
        if ~exist('fh','var') && warn==1
            warnfile = [filename(1:end-4),'-warnings.txt'];
            fh = fopen(warnfile,'wt');
      
        elseif warn==1            
            fprintf(fh, 'No particles found in frame %d\n',i);
            empty = []; % empty warn flag
        end
    else
        if intth > 0
            peaks = pkRefnd(im.filterImg(:,:,i),t_pkfind,boxr,intth);
            im.planeAttr(i).particleintthr = peaks;
        else
            peaks = t_pkfind;
            im.planeAttr(i).particleintthr = t_pkfind;
        end
        
        % Find subpixel resolution positions
        % (CJI 2008.08.10) recommend box size for cntrd is bpass
        % lengthscale+2 using recommended length scale could be an option
        t_cen = centfind(im.filterImg(:,:,i),peaks(:,1:2),boxr,bpsz,opt);
        % (BPB 2010.07.17) Delaunay requires doubles and it appears as if centfind returns
        % the centers in whatever class the input image is.
        t_cen = double(t_cen);
        t_npeaks = size(t_cen,1);
        
        % some outputs
        im.planeAttr(i).particlepos = t_cen;
        im.planeAttr(i).nparticle = t_npeaks;
        im.planeAttr(i).particlepxl = t_pkfind;
        
        if isempty(t_cen)
            
            if ~exist('fh','var') && warn==1
                warnfile = [filename(1:end-4),'-warnings.txt'];
                fh = fopen(warnfile,'wt');
                
            elseif warn==1
                fprintf(fh, 'No particles found in frame %d\n',i);
                empty = []; % empty warn flag
            end
            
            continue
        end
        % build for track.m
        t_pos = [t_pos;[t_cen(:,1:2),ones(t_npeaks,1)*i]];
                          
        % calculate the average interparticle distance this could be used
        % to set max_disp per frame.  Eliminating sptpara.max_disp
        if t_npeaks < 3
            
            % warn file
            if ~exist('fh','var') && warn==1
                warnfile = [filename(1:end-4),'-warnings.txt'];
                fh = fopen(warnfile,'wt');
                
            elseif warn==1
                fprintf(fh, 'Couldn''t calculate interparticle distance in frame %d\n',i);
                nodis = []; % nodis warn flag
            end
        else
            % Triangulation to find nearest neighbors
            TRI = delaunay(t_cen(:,1),t_cen(:,2));
            ntri = size(TRI,1);
            ii=1;
            for j = 1:ntri
                x= t_cen(TRI(j,:),1);
                dx = x - circshift(x,-1);
                y= t_cen(TRI(j,:),2);
                dy = y - circshift(y,-1);
                len(ii:ii+2)=sqrt(dx.^2+dy.^2);
                ii=ii+3;
            end
            ulen = unique(len);
            hlenbar = mean(ulen)/2;
            im.planeAttr(i).iparticledis(1) = mean(ulen);
            im.planeAttr(i).iparticledis(2) = std(ulen);
            
            % Warn if average distance is less than max disp skip this if
            % using spt_interactive you aren't going to use the warning
            % anyway
            if isfield(sptpara,'max_disp') && hlenbar<sptpara.max_disp
                if ~exist('fh','var') && warn==1
                    warnfile = [filename(1:end-4),'-warnings.txt'];
                    fh = fopen(warnfile,'wt');
                    
                elseif warn==1
                    fprintf(fh, ['Max displacement(%d)is greater than half',...
                        ' the mean interparticle distance(%5.2f) found in frame %d\n']...
                        ,[sptpara.max_disp, hlenbar,i]);
                    dist = []; % dist warn flag
                end              
            end           
        end
    end
end

% Return to previous warning state
warning(wstate1);

% Display Warnings
if exist('fh','var')
    if exist('empty','var')
        fprintf(1,'Warning: %s had frames without detected particles.\n',filename);
        disp(['See ',warnfile,' for frame numbers']);
        clear empty
    end
    
    if exist('dist','var')
        fprintf(1,'Warning: %s had frames with a high density of particles.\n',filename);
        disp(['See ',warnfile,' for frame numbers']);
        clear dist
    end
    
    if exist('nodis','var')
        fprintf(1,'Warning: Couldn''t calculated interparticle distance in some frames of %s.\n',filename);
        disp(['See ',warnfile,' for frame numbers']);
        clear nodis
    end
    
    fclose(fh);
    clear fh
end
poslist = t_pos;
