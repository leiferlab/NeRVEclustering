% spt_calc will calculate the trajectories of single particles in a stack
% tiff movie.  It is made to be used with spt.m but can be used on its on
% as well. Note: Some error checking is done in spt.m wrapper so only use
% this if you know what you are doing :-)

function [trajlist,traj] = spt_track(sptpara,poslist)

% INPUTS:
% sptpara.mtl = minimum trajectory length: a trajectory must at least be
%               this length in frames to be retained.
%
% sptpara.max_disp = maximum particle displacement:  this is used by
%                    track.m to set a maximum search region for particle
%                    linking.  Should be set to slightly less than the
%                    average interparticle distance but must be >> than RMS
%                    displacment between frames.  This is accomplished by
%                    keeping the particle density low.  A warning will be
%                    issued if 2*max_disp > average particle distance.
%                    Since the data begains to become questionable.  As
%                    max_disp approaches the average particle distance the
%                    amount of time needed to calcuate traj ballons and
%                    sometimes track will fail. 
%
% sptpara.trackMem = the number of frames that a particle can be "off" for
% to still come back as the "same" particle, passed into track as the "mem"
% parameter

% OUTPUTS:
% trajlist = an array contain list of traj found [xi, yi, idf, idt], xi
%            is is the xi, yi is the coordinate of the particle which is 
%            located in frame idf and belongs to trajectory idt.
%
% Dependencies: track.m
%
% Author: Colin Ingram
% Created: May 2010
%
% Revisions:
% Version     Date  Author   Description
%          2010.02.28 SB:  handles frames without any particles correctly
%% Inputs
mtl = sptpara.mtl;
max_disp = sptpara.max_disp;
t_pos = poslist;

%% Deal with non-consective frames


% detecting consequtive frames and definfing a new matrix to help the
% segregation of t_pos based on that
lastpos=find(t_pos(:,3),1,'last');
t_posnew=zeros(lastpos,1);
t_posnew(1)=0;
for i=1:lastpos-1
    if t_pos((i+1),3)==t_pos(i,3)% for frames with more than one particle
        t_posnew(i+1)=0;
    elseif t_pos((i+1),3)>t_pos(i,3)
        t_posnew(i+1)=t_pos((i+1),3)-t_pos(i,3)-1;
    else
    end
end
%Track
% allow other code to let particles blink
if not(isfield(sptpara,'trackMem'))
    param.mem   = 0;
else
    param.mem = sptpara.trackMem;
end
param.dim   = 2;
param.good  = mtl+1;
param.quiet = 1;
jump = find(t_posnew);
trajlist=[];
% track if all frames are consecitive
if isempty(jump)
    trajlist = track(t_pos,max_disp,param);
    
% break up tracklist and handoff seperatly
else
    sublenstemp = [(jump(1)-1),(diff(jump))'];
    sublens=[sublenstemp,(lastpos-max(jump)+1)];
    submats = mat2cell(t_pos, (sublens)', size(t_pos,2));
    nummat=size(submats); % numer of segments of t_pos
    count=0;
    for i=1:nummat(1)
        block=cell2mat(submats(i));
        a=max(block(:,3))- min(block(:,3));
        if a >=mtl+1 % selecting matrices of enough steps
%             mtl=3;
            try
                Trajlist=track(block,max_disp,param);
            catch ME
                fprintf(1,['less than ',num2str(mtl),' step traj. \ngoing to the next traj \n'])
                continue
            end
            for j=1:size(Trajlist,1)
                Trajlist(j,4)=Trajlist(j,4)+count;
            end
            count=max(Trajlist(:,4));
            trajlist=vertcat(trajlist,Trajlist);
        end
    end
    % modifying trajlist for the trajno
    if isempty(trajlist);
        warning('spt_track:noTrajectories','No trajectories could be made')
        traj = struct('pos',[],'length',0,'frame',0);
        return
    end
    trajcount=trajlist(:,4);
    for i=2:size(trajlist,1)
        if trajcount(i)-trajcount(i-1)>0
            trajlist(i,4)=trajlist((i-1),4)+1;
        else trajlist(i,4)=trajlist((i-1),4);
        end
    end
    trajlist(:,4)=trajlist(:,4)-min(trajlist(:,4))+1;
end

%traj = [];
%Renee added to avoid having to use spt_dispmoments
n_traj = max(trajlist(:,4));
traj(n_traj) = struct('pos',[],'length',0,'frame',0);

for i = 1:n_traj
    
     % Grab trajectory and set up variables
    t_slice = trajlist(trajlist(:,4) == i,:);
    t_nsteps = size(t_slice,1)-1;
    t_frame1 = t_slice(1,3);
    
     % Output
    traj(i).pos = t_slice(:,1:3);
    traj(i).length = t_nsteps;
    traj(i).frame = t_frame1;
    
end