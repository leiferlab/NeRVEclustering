function tform=makeTPStform(imageSize,model,moving)


    d=size(moving,2);
    subSamp=[50,50,5];
    [X,Y,Z]=ndgrid(1:subSamp(1):imageSize(1),1:subSamp(2):imageSize(2),1:subSamp(3):imageSize(3));
    
    for current_level = 1:level
       current_config = initialize_config(model, scene, 'tps', scales(current_level), lambdas(current_level), [], fiducial_indices, spring_constant,showFlag);
%%
scales=[5 .1 .1];lambdas=[0 0 0];
for current_level=1:3
               current_config = initialize_config(model, moving, 'tps', scales(current_level), lambdas(current_level), [], 1:8, 0,1);

        [current_param, model, history, config] = gmmreg_L2(current_config);
end

    %%    
       x=  transform_pointset([X(:) Y(:) Z(:)], 'tps', current_param, current_config.ctrl_pts,config.init_affine)
        x=reshape(x,size(X,1),size(X,2),size(X,3),3);
        Xout=x(:,:,:,1);Yout=x(:,:,:,2);Zout=x(:,:,:,3);
         
        [warped_pts, bending_energy] = transform_by_tps(current_param, current_config.ctrl_pts, current_config.ctrl_pts)
        if current_config.opt_affine == 0
            current_param = [repmat([zeros(1,d) 1],1,d), current_param];
        end
        multilevel_ctrl_pts = cat(3,  multilevel_ctrl_pts, current_config.ctrl_pts);
        multilevel_param = cat(3, multilevel_param, current_param);
    end













[K,U] = compute_kernel(ctrl_pts,moving);
    moving1=[moving ones(size(moving(:,1)))];
KPmatrix=[K moving1 ; moving1' zeros(size(moving1,2))];
moving2=[moving1 ;zeros(size(moving1,2))]
ctrl_pnts2=0
Wparams=0

tform=[];
