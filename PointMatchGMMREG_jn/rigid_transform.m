%% this function performs a rigid transform (translation, rotation, optional scaling) on M that best matchs a model M to a scene S.
function [rigid_transformed_M] = rigid_transform(M, S)
    M_transpose = transpose(M);
    rigid_transformed_M = M; %make it the same size as M for speed
    answer = absor(M_transpose,transpose(S));
    transformation_matrix = getfield(answer,'M');
    for i = 1:size(rigid_transformed_M, 1)
        newCoord = transformation_matrix*[M_transpose(:, i); 1];
        newCoord = transpose(newCoord);
        newCoord(:,4) = [];
        rigid_transformed_M(i, :) = newCoord();
    end
end
