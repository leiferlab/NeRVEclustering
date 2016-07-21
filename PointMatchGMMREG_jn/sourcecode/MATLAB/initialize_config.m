function [config] = initialize_config(model, scene, motion, scale, lambda, init_param, fiducial_indices, spring_constant,showFlag)

if nargin==8;
    showFlag=0;
end


config.model = model;
config.scene = scene;
config.motion = motion;
% estimate the scale from the covariance matrix
[n,d] = size(model);
if d>3
    d=3;
end
    
config.scale = power(det(model'*model/n), 1/(2^d));
config.display = showFlag;
config.init_param = [ ];
config.max_iter = 400;
config.normalize = 0;
config.fiducial_indices = fiducial_indices;
config.spring_constant = spring_constant;
switch lower(motion)
    case 'tps'
        interval = 5;
        config.ctrl_pts =  set_ctrl_pts(model, scene, interval);
        config.alpha = 1;
        config.opt_affine = 0;
        [n,d] = size(config.ctrl_pts); % number of points in model set
        if nargin > 3
            %optional scale and lambda parameters are specified
            config.scale = scale;
            config.beta = lambda;
        else
            config.beta = 0.08; %same as lambda
        end
        if size(init_param,1) > 0
            %the initial parameter is specified (not [])
            p = reshape(init_param, d, size(init_param,2)/d); p = p';
            config.init_tps = p(end+1-d*(n-d-1):end);
            config.init_affine = p(1:d*(d+1));
            config.init_param = init_param;
        else
            config.init_tps = zeros(n-d-1,d);
            init_affine = repmat([zeros(1,d) 1],1,d);
            config.init_param = [init_affine zeros(1, d*n-d*(d+1))];
            config.init_affine = [ ];
        end
    otherwise
        [x0,Lb,Ub] = set_bounds(motion);
        config.init_param = x0;
        config.Lb = Lb;
        config.Ub = Ub;
end

