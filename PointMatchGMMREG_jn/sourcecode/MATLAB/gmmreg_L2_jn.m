%function [param, tt] = GMMReg(model, scene, scale, motion, display, init_param);
%   'model' and 'scene'  are two point sets
%   'scale' is a free scalar parameter
%   'motion':  the transformation model, can be
%         ['rigid2d', 'rigid3d', 'affine2d', 'affine3d']
%         The default motion model is 'rigid2d' or 'rigid3d' depending on
%         the input dimension
%   'display': display the intermediate steps or not.
%   'init_param':  initial parameter

function [param, transformed_model, history, config] = gmmreg_L2_jn(config)
%%=====================================================================
%% $Author: bing.jian $
%% $Date: 2009-02-10 02:13:49 -0500 (Tue, 10 Feb 2009) $
%% $Revision: 121 $
%%=====================================================================

% todo: use the statgetargs() in statistics toolbox to process parameter name/value pairs
% Set up shared variables with OUTFUN
history.x = [ ];
history.fval = [ ];
if nargin<1
    error('Usage: gmmreg_L2(config)');
end
[n,d] = size(config.model); % number of points in model set
if d>3
    d=3;
end

if (d~=2)&&(d~=3)
    error('The current program only deals with 2D or 3D point sets.');
end

options = optimset( 'display','off', 'LargeScale','off','GradObj','on', 'TolFun',1e-10, 'TolX',1e-010, 'TolCon', 1e-10);
options = optimset(options, 'outputfcn',@outfun);
options = optimset(options, 'MaxFunEvals', config.max_iter);
options = optimset(options, 'GradObj', 'on');

switch lower(config.motion)
    case 'tps'
        scene = config.scene;
        scale = config.scale;
        alpha = config.alpha;
        beta = config.beta;
        model=config.model;
        [n,d] = size(config.ctrl_pts);
        [m,d] = size(config.model);
        if d>3
            d=3;
        end
        
        [K,U] = compute_kernel(config.ctrl_pts, config.model(:,1:d));
        Pm = [ones(m,1) config.model(:,1:d)];
        Pn = [ones(n,1) config.ctrl_pts(:,1:d)];
        PP = null(Pn');  % or use qr(Pn)
        basis = [Pm U*PP];
        kernel = PP'*K*PP;

        init_tps = config.init_tps;  % it should always be of size d*(n-d-1)
        if isempty(config.init_affine)
            % for your convenience, [] implies default affine
            config.init_affine = repmat([zeros(1,d) 1],1,d);
        end
        if config.opt_affine % optimize both affine and tps
            init_affine = [ ];
            x0 = [config.init_affine init_tps(end+1-d*(n-d-1):end)];
        else % optimize tps only
            init_affine = config.init_affine;
            x0 = init_tps(end+1-d*(n-d-1):end);
        end
             %           [param,Emin,exitflag] = fminunc(@(x)gmmreg_L2_tps_costfunc(x, init_affine, basis, kernel, scene, scale, alpha, beta, n, d, config.fiducial_indices, config.spring_constant), x0,  options);
%if extra fiducial data is present, use it, otherwise, do not
             if size(model,2)>d
            [param,Emin,exitflag] = fminunc(@(x)gmmreg_L2_tps_costfunc_jn(x, init_affine, basis, kernel, scene, scale, alpha, beta, n, d, config.fiducial_indices, config.spring_constant,model), x0,  options);
else
            [param,Emin,exitflag] = fminunc(@(x)gmmreg_L2_tps_costfunc_jn(x, init_affine, basis, kernel, scene, scale, alpha, beta, n, d, config.fiducial_indices, config.spring_constant), x0,  options);

end

            % param2=param;
% for i=1:100
%             [E(i),Egrad]=gmmreg_L2_tps_costfunc_jn(param2, init_affine, basis, kernel, ...
%     scene, scale, alpha, beta, n, d, config.fiducial_indices,...
%     config.spring_constant);
% param2=param2-Egrad/10000;
%                            transformed = transform_pointset(config.model, config.motion, param2, config.ctrl_pts, init_affine);     
%                    dist = L2_distance(transformed,config.scene,config.scale);
%                    DisplayPoints(transformed,config.scene,d);
%                    drawnow
%                        title(sprintf('L2distance: %f',dist));
%                    drawnow;
% end

            
            
 
        transformed_model = transform_pointset(config.model(:,1:d), config.motion, param, config.ctrl_pts(:,1:d), init_affine);
if size(model,2)>d
    transformed_model=[transformed_model model(:,end-1:end)];
end
        
        if config.opt_affine
            config.init_tps = param(end+1-d*(n-d-1):end);
            config.init_affine = param(1:d*(d+1));
        else
            config.init_tps = param;
        end
    otherwise
        x0 = config.init_param;
        param = fmincon(@gmmreg_L2_costfunc, x0, [ ],[ ],[ ],[ ], config.Lb, config.Ub, [ ], options, config);
        transformed_model = transform_pointset(config.model, config.motion, param);
        config.init_param = param;
end

    function stop = outfun(x,optimValues,state,varargin)
     stop = false;
     if config.display
     displayCounter=getappdata(0,'displaycounter');
     if isempty(displayCounter)
         displayCounter=0;
     else
         displayCounter=displayCounter+1;
     end
          setappdata(0,'displaycounter',round(displayCounter));
     else
         displayCounter=0;
     end
     switch state
         case 'init'
             if config.display>0
               set(gca,'FontSize',16);
             end
         case 'iter'
            if config.display>0 && ~mod(displayCounter,config.display)

               history.fval = [history.fval; optimValues.fval];
               history.x = [history.x; reshape(x,1,length(x))];
                   hold off
                   switch lower(config.motion)
                       case 'tps'
                           transformed = transform_pointset(config.model(:,1:d), config.motion, x, config.ctrl_pts(:,1:d), init_affine);
                       otherwise
                           transformed = transform_pointset(config.model(:,1:d), config.motion, x);
                   end
                   dist = L2_distance(transformed,config.scene,config.scale);
                   DisplayPoints(transformed(:,1:d),config.scene(:,1:d),d);
                   title(sprintf('L2distance: %f',dist));
                   drawnow;
               end
         case 'done'
              %hold off
         otherwise
     end
    end


end



function [dist] = L2_distance(model, scene, scale)
    dist = gaussOverlapSelf(model,scale) + gaussOverlapSelf(scene,scale) - 2*gaussOverlap(model,scene,scale);
end


