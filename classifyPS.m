function track_matrix=classifyPS(varargin)
%classifyPS takes an input pointStats structure from a single volume and
%adds the trackIdx for each of the neurons and the track weights.

% Inputs:
% classifyPS(PS,dataFolder)
%   PS-pointStats file for a single volume, this is the output of the
%   WormStraightening codes.
%   dataFolder - dataFolder (normally BrainScanner) that contains a
%   PS_ref.mat file and PointsStats2_info.mat file.

% classifyPS(PS,PS_ref,masterVec);
%   PS_ref - the PS_ref file containg the pointstats data for the reference
%   volumes. Normally made form the makePointStatsRef.m program.
%   masterVec - the centers of each of the identified clusters, normally
%   from clusterWormTrackCompiler code, saved in the _info.mat file.

% classifyPS(PS,PS_ref,masterVec,hitCutoff)
%   hitCutoff - threshold used for definiing membership to cluster,

%parse inputs
if nargin==2
    dataFolder=varargin{2};
    PS_ref=load([dataFolder filesep 'pointStatsRef.mat']);
    PS_ref=PS_ref.PS_ref;
    PSref_info=load([dataFolder filesep 'PointsStats2_info.mat']);
    masterVec=PSref_info.masterVec;
    if isfield(PSref_info,'hitCutoff')
        hitCutoff=PSref_info.hitCutoff;
    else
        hitCutoff=0.2;
    end
elseif nargin>=3
    PS_ref=varargin{2};
    masterVec=varargin{3};
    if  nargin==4
        hitCutoff=varargin{4};
    else
        hitCutoff=0.2;
    end
endtrack_matrix=compareWithRef(PS,PS_ref);
PS=varargin{1};

%make track_matrix by doing the registration between the points in PS and
%the points in PS_ref

track_matrix=compareWithRef(PS,PS_ref);

% % count max number of neurons
% n_ref_neurons=cellfun(@(x) size(x,1),{PS_ref.straightPoints});
% % do one hot encoding to turn matches into binary matrix
% track_matrix_bin=oneHotNeuron(track_matrix,n_ref_neurons);
% 
% %classify points based on distance to centers in the masterVec
% output=distanceClassify(masterVec,track_matrix_bin,hitCutoff);
% 
% %save outputs back into pointStats structure
% PS.trackIdx=output{1};
% PS.trackWeights=output{2};


