function [energy, grad] = gmmreg_L2_tps_costfunc_jn(param, init_affine, basis,  kernel, scene, scale, alpha, beta, n, d, fiducial_indices, spring_constant,model)
    %%=====================================================================
    %% $RCSfile: gmmreg_L2_tps_costfunc.m,v $
    %% $Author: bing.jian $
    %% $Date: 2009-02-10 02:13:49 -0500 (Tue, 10 Feb 2009) $
    %% $Revision: 121 $
    %%=====================================================================
    alpha=alpha*10;
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

    if exist('model','var')
        after_tps=[after_tps model(:,end-1:end)];
    end
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
end
function [f, g] = general_costfunc(A, B, scale)

    [f1, g1] = gaussOverlap(A,A,scale);
    [f2, g2] = gaussOverlap(A,B,scale);
    f =  f1 - 2*f2;
    g = 2*g1 - 2*g2; 
end


function [f, g] = fiducial_costfunc(model, scene, fiducial_indices, scale, spring_constant) %pass in a list of which neurons are fiducials
%remove the fiducial neurons from the gausstransform
if any(fiducial_indices)
    model_excluding_fiducials=model(fiducial_indices(:,1),:);
    scene_excluding_fiducials =  scene(fiducial_indices(:,2),:);
else
end
%JN: just usiung all points for now
[f1, g1, f2, g2]=gaussOverlapDouble(model,scene,scale);
non_fiducial_g = 2*g1 - 2*g2; %calculate the non-fiducial gradient

g = non_fiducial_g(:,1:3); %the gradient array returned will be the same size as model, so initialize it as model
%calculate and reinsert the gradient for the fiducial neurons and also
    %calculate the sum of the spring energies to add to the final cost to
    %minimize
    if any(fiducial_indices)
    difference = model(fiducial_indices(:,1),:) - scene(fiducial_indices(:,2),:);
        g(:,1:3) = spring_constant * difference+non_fiducial_g;
        spring_energy = (0.5 * spring_constant * sum(difference(:).^2));
    else
        spring_energy=0;
    end
    f =  f1 - 2*f2+spring_energy; %the total energy that needs to be minimized is the gaussian transform energy plus the spring energy

end
