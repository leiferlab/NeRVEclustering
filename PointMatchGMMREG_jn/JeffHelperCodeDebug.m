
%% finding the best scene out of the 126 worms
zscaling = [1, 0, 0; 0, 1, 0; 0, 0, 9.45];
scores = [];
                        
tic

best_perfect_scene_index = 0;
best_perfect_scene_count = 0;
best_avg_score_index = 0;
best_avg_score = 0;
best_average_distance_index = 0;
best_average_distance = 1000;
%fiducial_indices = [1,1;3,3;5,5;7,7;9,9;11,11];
fiducial_indices=[];

spring_constant = 10;

%%

for scene_index = 24:24  
%for scene_index = 1:size(newWormMasterPre,3)
    perfect_scene_count = 0;
    total_score = 0;
    total_distance = 0;
    %for model_index = 1:size(newWormMasterPre,3)
    for model_index = 1:size(newWormMasterPre,3)
    %for model_index = 1:size(newWormMasterPre,3)
        if scene_index ~= model_index
            M = newWormMasterPre(:,:, model_index);
            S = newWormMasterPre(:,:, scene_index);

            %correct for z
            for i = 1:size(M,1)
                M(i,:) = transpose(zscaling * transpose(M(i,:)));
            end
            for i = 1:size(S,1)
                S(i,:) = transpose(zscaling * transpose(S(i,:)));
            end
%%
            [score_avgdist, match_result, Transformed_M, Transformed_S, multilevel_ctrl_pts, multilevel_param] = run_gmmreg(M, S, S, [], [], fiducial_indices, spring_constant);
            score = score_avgdist(1);
            average_distance = score_avgdist(2);


           
                %display the plot
                [n,d,l] = size(multilevel_ctrl_pts); 
                ctrl_pts = multilevel_ctrl_pts(:,:,1);
                transformed_ctrl_pts = ctrl_pts;
                %%
                for level = 1:l
                    %loop through each level and repeat the tps transformation
                    %on the control points and save them in
                    %transformed_ctrl_pts
                    param = multilevel_param(:,:,level);
                    p = reshape(param, 3, size(param,2)/3); p = p';
                    transformed_ctrl_pts = cat(3, transformed_ctrl_pts, transform_by_tps(p, transformed_ctrl_pts(:,:,1), multilevel_ctrl_pts(:,:,level)));
                end

                subplot(1,2,1);
                DisplayPoints(M, S, size(M,2), 0, determine_border(M, S)); axis off; %display the original
                subplot(1,2,2);
                %display the transformed
                DisplayPoints(Transformed_M, Transformed_S, size(Transformed_S,2),  0, determine_border(Transformed_M, Transformed_S), transformed_ctrl_pts, match_result); axis off 
                [scene_index, model_index]
                score_avgdist
                pause
            
            total_score = total_score + score;
            total_distance = total_distance + (average_distance*score);    
        end
    end

    if perfect_scene_count > best_perfect_scene_count
        best_perfect_scene_count = perfect_scene_count;
        best_perfect_scene_index = scene_index;
    end

    if total_score / (size(newWormMasterPre,3)-1) > best_avg_score
        best_avg_score = total_score / (size(newWormMasterPre,3)-1);
        best_avg_score_index = scene_index;
    end

    if total_distance / total_score  < best_average_distance
        best_average_distance = total_distance / total_score;
        best_average_distance_index = scene_index;
    end
    scores = [scores; perfect_scene_count, total_score / (size(newWormMasterPre,3)-1), total_distance / total_score]
    scene_index
end

toc
best_perfect_scene_index
best_perfect_scene_count
best_avg_score_index
best_avg_score
best_average_distance_index
best_average_distance
