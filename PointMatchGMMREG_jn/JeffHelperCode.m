
%% takes 2 dimensional array of positions with neuron indicies (see array names) to create 3 dimensional array with the 3rd dimension being different worms
concatenatedMarkerFile3D = [];
for i = 1:(size(concatenatedmarkerfileFixed, 1)/12)
    wormOrd2 = zeros(12,3);
    for j = (i*12-11):(i*12)
        wormOrd2(-concatenatedmarkerfileFixed(j,4), :) = concatenatedmarkerfileFixed(j,1:3);
    end
    concatenatedMarkerFile3D = cat(3,concatenatedMarkerFile3D, wormOrd2);
end

%% ML code for rigid transformation of raw files
for i = 1:size(newWormMasterPre,3)
    currentWorm = transpose(concatenatedMarkerFile3D(:,:,i));
    answer = absor(currentWorm,transpose(newWormMaster(:,:,1)));
    M = getfield(answer,'M');
    for j = 1:size(currentWorm, 2)
        newCoord = M*[currentWorm(:, j); 1];
        newCoord = transpose(newCoord);
        newCoord(:,4) = [];
        concatenatedMarkerFile3DTransformed(j, :, i) = newCoord();
    end
end

%newWormMaster = cat(3,wormMaster, concatenatedMarkerFile3DTransformed);

%% creates the txt files for use by gmmreg .ini files
for i = 1:size(newWormMaster,3)
    dlmwrite(strcat('.\odr2_data\worm', num2str(i), '.txt'),newWormMaster(:,:,i), ' ');
end

%% getting all the worms that needs to be plotted
plotWorms = [];
for i = 1:size(newWormMaster,3)
    plotWorms = [plotWorms ; newWormMaster(:,:,i)];
end

%% % getting the color and size arrays for scatter plot
plotColorArray = [];
for i = 1:size(newWormMaster,3)
    plotColorArray = [plotColorArray ; colorarray];
end

plotSizeArray = [];
for i = 1:size(newWormMaster,3)
    plotSizeArray = [plotSizeArray , sizearray];
end
%% plotting scatter
%scatter3(plotWorms(:,1),plotWorms(:,2),plotWorms(:,3),plotSizeArray,plotColorArray,'fill')

%% plotting each neuron in scatter
plotWorms = [];
for i = 1:size(newWormMaster,3)
    plotWorms = [plotWorms ; newWormMaster(12,:,i)];
end

scatter3(plotWorms(:,1),plotWorms(:,2),plotWorms(:,3),300,'r','fill')

%% loading the paneuronal data
i = 1;
input_file_name = strcat('\\leiferdata.princeton.edu\data\PanNeuronal\20140508\his58tdtomato6\raw data\stackData\stack', sprintf('%04d',i), 'data.mat');

while exist(input_file_name, 'file')
    load(input_file_name, 'centroids')
    output_file_name = strcat('.\panneuronal\stack', sprintf('%04d',i), '.txt');
    dlmwrite(output_file_name, centroids, ' ');
    i = i+1;
    input_file_name = strcat('\\leiferdata.princeton.edu\data\PanNeuronal\20140508\his58tdtomato6\raw data\stackData\stack', sprintf('%04d',i), 'data.mat');
end
%% paneuronal
%methodarray = cellstr(['EM_TPS '; 'EM_GRBF'; 'TPS_L2 '; 'GRBF_L2'; 'TPS_KC '; 'GRBF_KC'; 'rigid  ']);
methodarray = cellstr(['rigid']);

for method_index = 1:size(methodarray,1)
    
    tic

    best_perfect_scene_index = 0;
    best_perfect_scene_count = 0;
    best_avg_score_index = 0;
    best_avg_score = 0;
    best_average_distance_index = 0;
    best_average_distance = 1000;
    
%     for scene_index = 35:35  
    for scene_index = 1:1
        perfect_scene_count = 0;
        total_score = 0;
        total_distance = 0;
%        for model_index = 1:size(newWormMaster,3)
        frame_track_output = trackOutput(trackOutput(:,8) == scene_index, :); %get the tracked output for this frame
        scene_key = zeros(max(frame_track_output(:,7)),1);
        for i = 1:size(frame_track_output,1)
            scene_key(frame_track_output(i,7)) = frame_track_output(i,9);
        end
        
        for model_index = 3:900
            if scene_index ~= model_index
                frame_track_output = trackOutput(trackOutput(:,8) == model_index, :); %get the tracked output for this frame
                model_key = zeros(max(frame_track_output(:,7)),1);
                for i = 1:size(frame_track_output,1)
                    model_key(frame_track_output(i,7)) = frame_track_output(i,9);
                end
                
                [score_avgdist, match_result] = run_gmmreg(model_index, scene_index, char(methodarray(method_index)), scene_key, model_key);
                score = score_avgdist(1);
                average_distance = score_avgdist(2);
                %display the plot
                model_file = ml_GetPrivateProfileString('FILES','model', '.\panneuronal.ini');
                scene_file = ml_GetPrivateProfileString('FILES','scene', '.\panneuronal.ini');
                ctrl_pts_file = ml_GetPrivateProfileString('FILES','ctrl_pts', '.\panneuronal.ini');
                transformed_file = ml_GetPrivateProfileString('FILES','transformed_model', '.\panneuronal.ini');
                final_affine_file = ml_GetPrivateProfileString('FILES','final_affine', '.\panneuronal.ini');
                final_tps_file = ml_GetPrivateProfileString('FILES','final_tps', '.\panneuronal.ini');

                M = load(model_file);
                S = load(scene_file);
                ctrl_pts = load(ctrl_pts_file);
                Transformed_M = load(transformed_file);
                final_tps = load(final_tps_file);
                final_affine = load(final_affine_file);


                [n,d] = size(ctrl_pts);
                p = cat(1, final_affine, final_tps);
                Transformed_ctrl_pts = transform_by_tps(p, ctrl_pts, ctrl_pts);

                subplot(1,2,1);
                DisplayPoints(M, S, size(M,2), 0, determine_border(M, S)); axis off;
                subplot(1,2,2);
                DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts - ctrl_pts); axis off
                %DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts); axis off
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

        if total_score / (size(newWormMaster,3)-1) > best_avg_score
            best_avg_score = total_score / (size(newWormMaster,3)-1);
            best_avg_score_index = scene_index;
        end

        if total_distance / total_score < best_average_distance
            best_average_distance = total_distance / total_score;
            best_average_distance_index = scene_index;
        end
        %scene_index
    end

    toc
    char(methodarray(method_index))
    best_perfect_scene_index
    best_perfect_scene_count
    best_avg_score_index
    best_avg_score
    best_average_distance_index
    best_average_distance
end


%% finding the best sigma (scale)
%methodarray = cellstr(['EM_TPS '; 'EM_GRBF'; 'TPS_L2 '; 'GRBF_L2'; 'TPS_KC '; 'GRBF_KC'; 'rigid  ']);
methodarray = cellstr(['TPS_L2']);

possible_sigma_values = [5, 1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005];
sigma_combinations = zeros(size(possible_sigma_values, 2)^3, 6);
sigma_combinations_index = 1;
for sigma_index_1 = 1:size(possible_sigma_values,2)
    for sigma_index_2 = 1:size(possible_sigma_values,2)
        for sigma_index_3 = 1:size(possible_sigma_values,2)
            sigma_combinations(sigma_combinations_index,:) = [possible_sigma_values(sigma_index_1) possible_sigma_values(sigma_index_2) possible_sigma_values(sigma_index_3) 0 0 0];
            sigma_combinations_index = sigma_combinations_index + 1;
        end
    end
end

for method_index = 1:size(methodarray,1)
    
    tic

    best_perfect_sigma_index = 0;
    best_perfect_sigma_count = 0;
    best_avg_score_index = 0;
    best_avg_score = 0;
    best_average_distance_index = 0;
    best_average_distance = 1000;
    scene_index = 35;
    
    
    for sigma_index = 1:size(sigma_combinations,1) 
%    for scene_index = 1:size(newWormMaster,3)
        perfect_scene_count = 0;
        total_score = 0;
        total_distance = 0;
%        for model_index = 1:size(newWormMaster,3)
        for model_index = 1:size(newWormMaster,3)
            if scene_index ~= model_index
                sigmastr = [num2str(sigma_combinations(sigma_index,1)), ' ', num2str(sigma_combinations(sigma_index,2)), ' ', num2str(sigma_combinations(sigma_index,3))];
                [score_avgdist, match_result] = run_gmmreg(model_index, scene_index, char(methodarray(method_index)), [], [], sigmastr);
                score = score_avgdist(1);
                average_distance = score_avgdist(2);
                if score == 12
                    %perfect
                    perfect_scene_count = perfect_scene_count + 1;
                else
%                     display the plot
%                     model_file = ml_GetPrivateProfileString('FILES','model', '.\odr2.ini');
%                     scene_file = ml_GetPrivateProfileString('FILES','scene', '.\odr2.ini');
%                     ctrl_pts_file = ml_GetPrivateProfileString('FILES','ctrl_pts', '.\odr2.ini');
%                     transformed_file = ml_GetPrivateProfileString('FILES','transformed_model', '.\odr2.ini');
%                     final_affine_file = ml_GetPrivateProfileString('FILES','final_affine', '.\odr2.ini');
%                     final_tps_file = ml_GetPrivateProfileString('FILES','final_tps', '.\odr2.ini');
% 
%                     M = load(model_file);
%                     S = load(scene_file);
%                     ctrl_pts = load(ctrl_pts_file);
%                     Transformed_M = load(transformed_file);
%                     final_tps = load(final_tps_file);
%                     final_affine = load(final_affine_file);
% 
% 
%                     [n,d] = size(ctrl_pts);
%                     p = cat(1, final_affine, final_tps);
%                     Transformed_ctrl_pts = transform_by_tps(p, ctrl_pts, ctrl_pts);
% 
%                     subplot(1,2,1);
%                     DisplayPoints(M, S, size(M,2), 0, determine_border(M, S)); axis off;
%                     subplot(1,2,2);
%                     DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts - ctrl_pts, match_result); axis off
%                     %DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts); axis off
%                     [scene_index, model_index]
%                     score_avgdist
%                     pause
                end
                total_score = total_score + score;
                total_distance = total_distance + (average_distance*score);
            end
        end
        
        sigma_combinations(sigma_index, 4) = perfect_scene_count;
        sigma_combinations(sigma_index, 5) = total_score / (size(newWormMaster,3)-1);
        sigma_combinations(sigma_index, 6) = total_distance / total_score;
        
        if perfect_scene_count > best_perfect_sigma_count
            best_perfect_sigma_count = perfect_scene_count;
            best_perfect_sigma_index = sigma_index;
        end

        if total_score / (size(newWormMaster,3)-1) > best_avg_score
            best_avg_score = total_score / (size(newWormMaster,3)-1);
            best_avg_score_index = sigma_index;
        end

        if total_distance / total_score < best_average_distance
            best_average_distance = total_distance / total_score;
            best_average_distance_index = sigma_index;
        end
        sigma_index
        perfect_scene_count
    end

    toc
    char(methodarray(method_index))
    sigma_combinations
    best_perfect_sigma_index
    best_perfect_sigma_count
    best_avg_score_index
    best_avg_score
    best_average_distance_index
    best_average_distance
end


%% finding the best lambda
%methodarray = cellstr(['EM_TPS '; 'EM_GRBF'; 'TPS_L2 '; 'GRBF_L2'; 'TPS_KC '; 'GRBF_KC'; 'rigid  ']);
methodarray = cellstr(['TPS_L2']);

possible_lambda_values = [0.08];
lambda_combinations = zeros(size(possible_lambda_values, 2)^3, 6);
lambda_combinations_index = 1;
for lambda_index_1 = 1:size(possible_lambda_values,2)
    for lambda_index_2 = 1:size(possible_lambda_values,2)
        for lambda_index_3 = 1:size(possible_lambda_values,2)
            lambda_combinations(lambda_combinations_index,:) = [possible_lambda_values(lambda_index_1) possible_lambda_values(lambda_index_2) possible_lambda_values(lambda_index_3) 0 0 0];
            lambda_combinations_index = lambda_combinations_index + 1;
        end
    end
end

for method_index = 1:size(methodarray,1)
    
    tic

    best_perfect_lambda_index = 0;
    best_perfect_lambda_count = 0;
    best_avg_score_index = 0;
    best_avg_score = 0;
    best_average_distance_index = 0;
    best_average_distance = 1000;
    scene_index = 35;
    
    
    for lambda_index = 1:size(lambda_combinations,1) 
%    for scene_index = 1:size(newWormMaster,3)
        perfect_scene_count = 0;
        total_score = 0;
        total_distance = 0;
%        for model_index = 1:size(newWormMaster,3)
        for model_index = 1:size(newWormMaster,3)
            if scene_index ~= model_index
                lambdastr = [num2str(lambda_combinations(lambda_index,1)), ' ', num2str(lambda_combinations(lambda_index,2)), ' ', num2str(lambda_combinations(lambda_index,3))];
                [score_avgdist, match_result] = run_gmmreg(model_index, scene_index, char(methodarray(method_index)), [], [], lambdastr);
                score = score_avgdist(1);
                average_distance = score_avgdist(2);
                if score == 12
                    %perfect
                    perfect_scene_count = perfect_scene_count + 1;
                else
%                     display the plot
%                     model_file = ml_GetPrivateProfileString('FILES','model', '.\odr2.ini');
%                     scene_file = ml_GetPrivateProfileString('FILES','scene', '.\odr2.ini');
%                     ctrl_pts_file = ml_GetPrivateProfileString('FILES','ctrl_pts', '.\odr2.ini');
%                     transformed_file = ml_GetPrivateProfileString('FILES','transformed_model', '.\odr2.ini');
%                     final_affine_file = ml_GetPrivateProfileString('FILES','final_affine', '.\odr2.ini');
%                     final_tps_file = ml_GetPrivateProfileString('FILES','final_tps', '.\odr2.ini');
% 
%                     M = load(model_file);
%                     S = load(scene_file);
%                     ctrl_pts = load(ctrl_pts_file);
%                     Transformed_M = load(transformed_file);
%                     final_tps = load(final_tps_file);
%                     final_affine = load(final_affine_file);
% 
% 
%                     [n,d] = size(ctrl_pts);
%                     p = cat(1, final_affine, final_tps);
%                     Transformed_ctrl_pts = transform_by_tps(p, ctrl_pts, ctrl_pts);
% 
%                     subplot(1,2,1);
%                     DisplayPoints(M, S, size(M,2), 0, determine_border(M, S)); axis off;
%                     subplot(1,2,2);
%                     DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts - ctrl_pts, match_result); axis off
%                     %DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts); axis off
%                     [scene_index, model_index]
%                     score_avgdist
%                     pause
                end
                total_score = total_score + score;
                total_distance = total_distance + (average_distance*score);
            end
        end
        
        lambda_combinations(lambda_index, 4) = perfect_scene_count;
        lambda_combinations(lambda_index, 5) = total_score / (size(newWormMaster,3)-1);
        lambda_combinations(lambda_index, 6) = total_distance / total_score;
        
        if perfect_scene_count > best_perfect_lambda_count
            best_perfect_lambda_count = perfect_scene_count;
            best_perfect_lambda_index = lambda_index;
        end

        if total_score / (size(newWormMaster,3)-1) > best_avg_score
            best_avg_score = total_score / (size(newWormMaster,3)-1);
            best_avg_score_index = lambda_index;
        end

        if total_distance / total_score < best_average_distance
            best_average_distance = total_distance / total_score;
            best_average_distance_index = lambda_index;
        end
        lambda_index
        perfect_scene_count
    end

    toc
    char(methodarray(method_index))
    lambda_combinations
    best_perfect_lambda_index
    best_perfect_lambda_count
    best_avg_score_index
    best_avg_score
    best_average_distance_index
    best_average_distance
end


%% finding the best scene out of the 126 worms
%methodarray = cellstr(['EM_TPS '; 'EM_GRBF'; 'TPS_L2 '; 'GRBF_L2'; 'TPS_KC '; 'GRBF_KC'; 'rigid  ']);
methodarray = cellstr(['TPS_L2 ']);
zscaling = [1, 0, 0; 0, 1, 0; 0, 0, 7];
scores = [];
                        
for method_index = 1:size(methodarray,1)
    
    tic

    best_perfect_scene_index = 0;
    best_perfect_scene_count = 0;
    best_avg_score_index = 0;
    best_avg_score = 0;
    best_average_distance_index = 0;
    best_average_distance = 1000;

%    for scene_index = 35:35  
    for scene_index = 1:size(newWormMasterPre,3)
        perfect_scene_count = 0;
        total_score = 0;
        total_distance = 0;
%        for model_index = 1:size(newWormMaster,3)
        for model_index = 1:size(newWormMasterPre,3)
            if scene_index ~= model_index
                M = newWormMasterPre(:,:, model_index);
                S = newWormMasterPre(:,:, scene_index);
                
                %correct for z
                if model_index > 35
                    for i = 1:size(M,1)
                        M(i,:) = transpose(zscaling * transpose(M(i,:)));
                    end
                end
                if scene_index > 35
                    for i = 1:size(S,1)
                        S(i,:) = transpose(zscaling * transpose(S(i,:)));
                    end
                end
                
                
                [score_avgdist, match_result] = run_gmmreg(M, S, rigid_scene_worm, char(methodarray(method_index)), [], []);
                score = score_avgdist(1);
                average_distance = score_avgdist(2);
                if score == 12
                    %perfect
                    perfect_scene_count = perfect_scene_count + 1;
                else
%                     %display the plot
%                     %model_file = ml_GetPrivateProfileString('FILES','model', '.\odr2_static.ini');
%                     scene_file = ml_GetPrivateProfileString('FILES','scene', '.\odr2_static.ini');
%                     ctrl_pts_file = ml_GetPrivateProfileString('FILES','ctrl_pts', '.\odr2_static.ini');
%                     transformed_file = ml_GetPrivateProfileString('FILES','transformed_model', '.\odr2_static.ini');
%                     final_affine_file = ml_GetPrivateProfileString('FILES','final_affine', '.\odr2_static.ini');
%                     final_tps_file = ml_GetPrivateProfileString('FILES','final_tps', '.\odr2_static.ini');
% 
%                     %M = load(model_file);
%                     Transformed_S = load(scene_file);
%                     ctrl_pts = load(ctrl_pts_file);
%                     Transformed_M = load(transformed_file);
%                     final_tps = load(final_tps_file);
%                     final_affine = load(final_affine_file);
% 
% 
%                     [n,d] = size(ctrl_pts);
%                     p = cat(1, final_affine, final_tps);
%                     Transformed_ctrl_pts = transform_by_tps(p, ctrl_pts, ctrl_pts);
% 
%                     subplot(1,2,1);
%                     DisplayPoints(M, S, size(M,2), 0, determine_border(M, S)); axis off;
%                     subplot(1,2,2);
%                     DisplayPoints(Transformed_M, Transformed_S, size(Transformed_S,2),  0, determine_border(Transformed_M, Transformed_S), ctrl_pts, Transformed_ctrl_pts - ctrl_pts, match_result); axis off
%                     %DisplayPoints(Transformed_M, S, size(S,2),  0, determine_border(M, S), ctrl_pts, Transformed_ctrl_pts); axis off
%                     [scene_index, model_index]
%                     score_avgdist
%                     pause
                end
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

        if total_distance / total_score < best_average_distance
            best_average_distance = total_distance / total_score;
            best_average_distance_index = scene_index;
        end
        scores = [scores; perfect_scene_count, total_score / (size(newWormMasterPre,3)-1), total_distance / total_score]
        scene_index
    end

    toc
    char(methodarray(method_index))
    best_perfect_scene_index
    best_perfect_scene_count
    best_avg_score_index
    best_avg_score
    best_average_distance_index
    best_average_distance
end
