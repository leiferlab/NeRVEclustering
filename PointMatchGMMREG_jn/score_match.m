%% returns the number of correct matches and the average distance of the correct matches in an array [score, avg_dist]
function [score_avg, match_result] = score_match(model, scene, threshold, model_key, scene_key)
    %match_result = set_match(model, scene, threshold);
    match_result = tracking_set_match(model, scene, threshold);
    totaldist = 0;
    score = 0;
    max_score = 0;
    if size(model_key, 1) > 0 && size(scene_key, 1) > 0
        %scene and model keys are provided
        for i = 1:size(match_result,1)
            if match_result(i,1) > 0 
                if scene_key(i) == model_key(match_result(i,1)) && scene_key(i) > 0 
                    totaldist = totaldist + match_result(i,2);
                    score = score + 1;
                end
                if scene_key(i) > 0 && model_key(match_result(i,1)) > 0
                    max_score = max_score + 1;
                end
            end
        end
    else
        for i = 1:12
        %for i = 1:size(match_result,1)
            if i == match_result(i,1)
                totaldist = totaldist + match_result(i,2);
                score = score + 1;
            end
            max_score = max_score + 1;
        end
    end

    if score > 0
        avg_dist = totaldist / score;
    else
        avg_dist = -1;
    end
    score_avg = [score, avg_dist, score / max_score];
end


%%returns a matrix match_result with the first column being the indecies of
%%points in the model that matches to the row numbers of the scene. The
%%second column is the distance between model and the scene point. This
%%method is no longer used
function [match_result] = set_match(model, scene,threshold)
    match_result = zeros(size(scene,1), 2);
    scene_index = 1;
    while (scene_index <= size(scene,1))
        %min_dist = norm(scene(scene_index,:) - model(1,:));
        if match_result(scene_index, 1) == 0
            %this scene has no matching point found for it yet
            distances = zeros(size(model,1), 1);
            for model_index = 1:size(model,1)
                %loop through all the points in model and find the closest one to
                %the point in scene
                distances(model_index) = norm(scene(scene_index,:) - model(model_index,:));
    %             if dist <= min_dist
    %                 % the euclidean distance is the lowest found so far, update
    %                 min_dist = dist;
    %                 min_dist_index = model_index;
    %             end
            end
            [sorted_distances, sorted_index] = sort(distances, 1, 'ascend'); %sort all the distances

            %  the model closest to the scene point should be sorted_index(1)
            for min_dist_index = 1:size(model,1)
                %loop though each of the minimum distances and assign scene to
                %the appropriate one
                if sorted_distances(min_dist_index) <= threshold
                    if ismember(sorted_index(min_dist_index), match_result(:,1))
                       %this model point has registered before with a scene point, register it with the
                       %closer one
                       previous_scene_index = find(match_result(:,1) == sorted_index(min_dist_index));
                       previous_dist = match_result(previous_scene_index, 2);
                       if previous_dist > sorted_distances(min_dist_index)
                           %newly found dist is smaller than previously found dist
                           match_result(scene_index, 1) = sorted_index(min_dist_index);
                           match_result(scene_index, 2) = sorted_distances(min_dist_index);

                           match_result(previous_scene_index, 1) = 0; %reset the previous found neuron as not found
                           match_result(previous_scene_index, 2) = 0;
                           scene_index = previous_scene_index - 1; %re-do from that scene again
                           break;
                       else
                           %newly found dist is larger than previously found
                           %dist, go to the next smallest distance
                       end
                    else
                       %this point has not registered before
                       match_result(scene_index, 1) = sorted_index(min_dist_index);
                       match_result(scene_index, 2) = sorted_distances(min_dist_index);
                       break;
                    end
                else
                    break; %there are no longer any values smaller than threshold, do not assign any model point to scene
                end
            end
        end

        scene_index = scene_index + 1;
    end
end


%%
function [match_result] = tracking_set_match(model, scene,threshold)
    globalmin = min(min(model(:)), min(scene(:))); %find the lowest number in both model and scene
    if (globalmin < 0)
        %there are negative numbers, so add it to the model and scene so
        %that tracking will work
        model = model - globalmin + 1;
        scene = scene - globalmin + 1;
    end

    scene_timepoint = cat(2, scene, zeros(size(scene,1),2)); %make scene time 0
    for i = 1:size(scene_timepoint,1)
        scene_timepoint(i, 4) = i; %make column 4 the index
    end
    model_timepoint = cat(2, model, ones(size(model,1),2)); %make model time 1
    for i = 1:size(model_timepoint,1)
        model_timepoint(i,4) = i; %make column 4 the index
    end
    all_points = cat(1, scene_timepoint, model_timepoint);

    track_match_result = track(all_points, threshold, struct('dim', 3, 'quiet', 1, 'difficult', 5.e+8, 'excessive', 2));
    match_result = zeros(size(scene,1),2);
    for i = 1:size(track_match_result, 1)
        if track_match_result(i, 5) == 0
            %the result is for a scene point (time 0)
           match_index = find(track_match_result(:,5) == 1 & track_match_result(:,6) == track_match_result(i,6));
           if match_index > 0
               %match found
               %track_match_result(i,4) is the scene index
               match_result(track_match_result(i,4),1) = track_match_result(match_index,4); %set the index of the matched point in model
               match_result(track_match_result(i,4),2) = norm(scene(track_match_result(i,4),:) - model(track_match_result(match_index,4),:)); %find the euclidean distance for the match
           end
        end
    end
end
