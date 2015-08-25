%%  gmmreg_L2_multilevel implements gmmreg algorithm going through several length scales
function [modelOut, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel_jn(model, scene, level, scales, lambdas, fiducial_indices, spring_constant,showFlag)
  % warps first input set to seconds, currently STRINCTLY 3D, this can be
  % changed though
  


multilevel_param = [];
    multilevel_ctrl_pts = [];
    if nargin==6;
        showFlag=0;
    end
    d=size(model,2);
    for current_level = 1:level
        tic
        current_config = initialize_config(model, scene, 'tps', scales(current_level), lambdas(current_level), [], fiducial_indices, spring_constant(current_level),showFlag);
        if current_level>1
            current_config.init_tps=config_out.init_tps;
        end
        
        [current_param, modelOut, ~, config_out] = gmmreg_L2_jn(current_config);
        if current_config.opt_affine == 0
            current_param = [repmat([zeros(1,d) 1],1,d), current_param];
        end
        multilevel_ctrl_pts = cat(3,  multilevel_ctrl_pts, current_config.ctrl_pts);
        multilevel_param = cat(3, multilevel_param, current_param);
        if showFlag
        display(['Finished ' num2str(current_level) ' in ' num2str(toc) 's']);
        end
    end
