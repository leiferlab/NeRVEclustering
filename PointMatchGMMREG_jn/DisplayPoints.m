function [axis_limits] = DisplayPoints(Model, Scene, dim, sampling, axis_limits, transformed_ctrl_pts, match_result)
%%=====================================================================
%% $RCSfile: DisplayPoints.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================

%if (nargin<4)
%    axis_limits = determine_border(Model,Scene);
%end

if (nargin<4)
    sampling = 0;
end

if (nargin<5)
    axis_limits = determine_border(Model, Scene);
end

if dim==2
    DisplayPoints2D(Model(:,1:2), Scene(:,1:2), sampling, axis_limits);
end

if dim==3
    if nargin > 5
        DisplayPoints3D(Model(:,1:3), Scene(:,1:3), sampling, axis_limits, transformed_ctrl_pts, match_result);
    else
        DisplayPoints3D(Model(:,1:3), Scene(:,1:3), sampling, axis_limits);
    end
end