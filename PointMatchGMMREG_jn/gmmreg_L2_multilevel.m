%%  gmmreg_L2_multilevel implements gmmreg algorithm going through several length scales
function [model, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(model, scene, level, scales, lambdas, fiducial_indices, spring_constant,showFlag)
  % warps first input set to seconds
  


multilevel_param = [];
    multilevel_ctrl_pts = [];
    if nargin==6;
        showFlag=0;
    end
    d=size(model,2);
    for current_level = 1:level
        current_config = initialize_config(model, scene, 'tps', scales(current_level), lambdas(current_level), [], fiducial_indices, spring_constant,showFlag);
        [current_param, model, ~, ~] = gmmreg_L2(current_config);
        if current_config.opt_affine == 0
            current_param = [repmat([zeros(1,d) 1],1,d), current_param];
        end
        multilevel_ctrl_pts = cat(3,  multilevel_ctrl_pts, current_config.ctrl_pts);
        multilevel_param = cat(3, multilevel_param, current_param);
    end
end