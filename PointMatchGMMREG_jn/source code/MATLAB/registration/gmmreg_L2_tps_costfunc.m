function [energy, grad] = gmmreg_L2_tps_costfunc(param, init_affine, basis,  kernel, scene, scale, alpha, beta, n, d, fiducial_indices, spring_constant)
    %%=====================================================================
    %% $RCSfile: gmmreg_L2_tps_costfunc.m,v $
    %% $Author: bing.jian $
    %% $Date: 2009-02-10 02:13:49 -0500 (Tue, 10 Feb 2009) $
    %% $Revision: 121 $
    %%=====================================================================

    if isempty(init_affine)
        %% if init_affine is given as [ ], then it means the affine matrix is 
        %% part of parameter and will be updated during optimization as well.
        %% In this case, the length of parameter should be n*d
        affine_param = reshape(param(1:d*(d+1)),d,d+1);
        affine_param = affine_param';
        tps_param = reshape(param(d*(d+1)+1:d*n),d,n-d-1);
        tps_param = tps_param';
    else
        %% if a non-empty init_affine is given, then it will be treated as
        %% a fixed affine matrix.
        %% In this case, the length of parameter should be (n-d-1)*d
        tps_param = reshape(param(1:d*n-d*(d+1)),d,n-d-1);
        tps_param = tps_param';
        affine_param = reshape(init_affine,d,d+1);
        affine_param = affine_param';
    end
    after_tps = basis*[affine_param;tps_param];
    bending = trace(tps_param'*kernel*tps_param);
    %[energy,grad] = general_costfunc(after_tps, scene, scale);
    [energy,grad] = fiducial_costfunc(after_tps, scene, fiducial_indices, scale, spring_constant);
    energy = alpha*energy + beta * bending;
    grad = alpha*basis'*grad;
    grad(d+2:n,:) = grad(d+2:n,:) + 2*beta*kernel*tps_param;
    if isempty(init_affine) 
        %% In this case, the length of gradient should be n*d    
        grad = grad';
        grad = reshape(grad,1,d*n);
    else 
        %% In this case, the length of parameter should be (n-d-1)*d    
        grad(1:d+1,:) = [ ];
        grad = grad';
        grad = reshape(grad,1,d*(n-d-1));
    end
    energy=energy
end
function [f, g] = general_costfunc(A, B, scale)

    [f1, g1] = GaussTransform(A,A,scale);
    [f2, g2] = GaussTransform(A,B,scale);
    f =  f1 - 2*f2;
    g = 2*g1 - 2*g2; 
end


function [f, g] = fiducial_costfunc(model, scene, fiducial_indices, scale, spring_constant) %pass in a list of which neurons are fiducials
%remove the fiducial neurons from the gausstransform
model_excluding_fiducials = zeros(size(model,1) - size(fiducial_indices,1), size(model,2));
model_excluding_fidicials_index = 1;
for i = 1:size(model, 1)
    if ~ismember(i, fiducial_indices(:,1))
        model_excluding_fiducials(model_excluding_fidicials_index,:) =  model(i,:);
    end
    model_excluding_fidicials_index = model_excluding_fidicials_index + 1;
end

scene_excluding_fiducials = zeros(size(scene,1) - size(fiducial_indices,1), size(model,2));
scene_excluding_fidicials_index = 1;
for i = 1:size(scene, 1)
    if ~ismember(i, fiducial_indices(:,2))
        scene_excluding_fiducials(scene_excluding_fidicials_index,:) =  scene(i,:);
    end
    scene_excluding_fidicials_index = scene_excluding_fidicials_index + 1;
end

%gauss transform the non-fiducial points
[f1, g1] = GaussTransform(model_excluding_fiducials,model_excluding_fiducials,scale);
[f2, g2] = GaussTransform(model_excluding_fiducials,scene_excluding_fiducials,scale);

non_fiducial_g = 2*g1 - 2*g2; %calculate the non-fiducial gradient

g = model; %the gradient array returned will be the same size as model, so initialize it as model
non_fiducial_g_index = 1;
spring_energy = 0;
%calculate and reinsert the gradient for the fiducial neurons and also
    %calculate the sum of the spring energies to add to the final cost to
    %minimize
for i = 1:size(g, 1)
    if ismember(i, fiducial_indices(:,1))
        %for the fiducial neurons,  their gradient is calculated as the derivative of the spring constant * the displacement
            %g(fiducial neuron)=k*[delta x, delta y, delt z]
        fiducial_index = find(fiducial_indices(:,1) == i); %get the index of the fiducial array
        difference = model(fiducial_indices(fiducial_index,1),:) - scene(fiducial_indices(fiducial_index,2),:);
        g(i,:) = spring_constant * difference;
        
        %for the fiducial neurons, add the cost term here with our energy function of model(fiducial) & scene(fiducial)
            %energy of a spring = 0.5 * k * x^2
            %+ 0.5 * k *(delta x ^2 + delta y^2 + delta z^2)
        spring_energy = spring_energy + (0.5 * spring_constant * sum(difference.^2));
    else
        %the non_fiducial gradient was already calculated
        g(i,:) = non_fiducial_g(non_fiducial_g_index, :);
        non_fiducial_g_index = non_fiducial_g_index + 1;
    end

end

f =  f1 - 2*f2 + spring_energy; %the total energy that needs to be minimized is the gaussian transform energy plus the spring energy

end
