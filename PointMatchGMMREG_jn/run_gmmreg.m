
%% runs the gmmreg for worms using various methods
function [score_avg, match_result, Transformed_M, S_rigid_transform, multilevel_ctrl_pts, multilevel_param] = run_gmmreg(M, S, rigid_transform_scene, model_key, scene_key, fiducial_indices, spring_constant)
    %use rigid transformation to align the two point sets
    M_rigid_transform = rigid_transform(M, rigid_transform_scene);
    S_rigid_transform = rigid_transform(S, rigid_transform_scene);
    
    %peform non-rigid transformation using TPS_L2 method
    [Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(M_rigid_transform, S_rigid_transform, 3, [1, 0.1, 0.01], [0.08, 0.08, 0.08], fiducial_indices, spring_constant);
    
    %socre the matches
    [score_avg, match_result] = score_match(Transformed_M, S_rigid_transform, 0.38, model_key, scene_key);
end
