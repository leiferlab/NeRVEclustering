function xyzs_avg = ActiveRevolutionFit(hyper_stack, cline_para, cline_initial)


% cline_para.alpha =10; 
% cline_para.beta = 200;

[row col, stack_z_size] = size(hyper_stack);

% initialize snake
xyzs = cline_initial;

%s_interp = 0:0.1:s(end);
%xyzs = interp1(s, cline_initial, s_interp, 'spline');
[m n] = size(xyzs);
    
% populating the penta diagonal matrix
A = zeros(m,m);
brow = zeros(1,m);
brow(1,1:5) = [cline_para.CLbeta, -(cline_para.CLalpha + 4*cline_para.CLbeta), (2*cline_para.CLalpha + 6 *cline_para.CLbeta), -(cline_para.CLalpha + 4*cline_para.CLbeta), cline_para.CLbeta];
A(1, 1:3) = [cline_para.CLbeta, -2*cline_para.CLbeta, cline_para.CLbeta];
A(2, 1:4) = [-cline_para.CLalpha-2*cline_para.CLbeta, 2*cline_para.CLalpha+5*cline_para.CLbeta, -cline_para.CLalpha-4*cline_para.CLbeta, cline_para.CLbeta];
for i=3:m-2
    A(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end
A(m-1, m-3:m) = [cline_para.CLbeta, -cline_para.CLalpha-4*cline_para.CLbeta, 2*cline_para.CLalpha+5*cline_para.CLbeta, -cline_para.CLalpha-2*cline_para.CLbeta];
A(m, m-2:m) = [cline_para.CLbeta, -2*cline_para.CLbeta, cline_para.CLbeta];
[L U] = lu(A + cline_para.gamma .* eye(m,m));
Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization

% find the center line of the cell



%%
Ar = zeros(m,m);
brow = zeros(1,m);
brow(1,1:5) = [cline_para.CLbetaR, -(cline_para.CLalphaR + 4*cline_para.CLbetaR), (2*cline_para.CLalphaR + 6 *cline_para.CLbetaR), -(cline_para.CLalphaR + 4*cline_para.CLbetaR), cline_para.CLbetaR];
Ar(1, 1:3) = [cline_para.CLbetaR, -2*cline_para.CLbetaR, cline_para.CLbetaR];
Ar(2, 1:4) = [-cline_para.CLalphaR-2*cline_para.CLbetaR, 2*cline_para.CLalphaR+5*cline_para.CLbetaR, -cline_para.CLalphaR-4*cline_para.CLbetaR, cline_para.CLbetaR];
for i=3:m-2
    Ar(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end
Ar(m-1, m-3:m) = [cline_para.CLbeta, -cline_para.CLalpha-4*cline_para.CLbeta, 2*cline_para.CLalpha+5*cline_para.CLbeta, -cline_para.CLalpha-2*cline_para.CLbeta];
Ar(m, m-2:m) = [cline_para.CLbeta, -2*cline_para.CLbeta, cline_para.CLbeta];
[L U] = lu(Ar + cline_para.gamma .* eye(m,m));
Arinv = inv(U) * inv(L); % Computing Ainv using LU factorization

% find the center line of the cell

%%


im = zeros(row, col, stack_z_size);
% two pass algorithm

% vidObj = VideoWriter('tmp.avi', 'Motion Jpeg AVI');
% vidObj.FrameRate = 15;
% vidObj.Quality = 95;
% open(vidObj);
    %k: stack_z_size index
        %obtain stack_z_size(t)

        %% filter image to "normalize" intensity
im=hyper_stack;
        %gradient
        if ismatrix(im)
            
        [Fy, Fx] = gradient(im*cline_para.gradient_force); %computing the gradient
       Fz=Fx*0;
        else
        
        [Fx, Fy, Fz] = gradient(im*cline_para.gradient_force); %computing the gradient
        end
        %% prepare for cline_para.iterations
        xyzs_avg = zeros(size(xyzs));
dtheta=pi/40;
theta=dtheta:dtheta:2*pi;
        counter = 0;
        %moving the snake in each iteration
        for i=1:cline_para.iterations;
            %calculate image gradient force on the nodes
                xs=xyzs(:,1);
                ys=xyzs(:,2);
                zs=xyzs(:,3);
                radii=abs(xyzs(:,4));
                radii=radii*0+26;
                [t_vector,b_vector,n_vector]=tbnVector([xs,ys,zs]);
b_vector=[b_vector(1,:);b_vector;];
n_vector=[n_vector(1,:);n_vector;];
signVector=sign(n_vector(:,1));
n_vector=bsxfun(@times,n_vector,signVector);


nx=(b_vector(:,1))*cos(theta)+(n_vector(:,1))*sin(theta);
ny=(b_vector(:,2))*cos(theta)+(n_vector(:,2))*sin(theta);
nz=(b_vector(:,3))*cos(theta)+(n_vector(:,3))*sin(theta);
rx=bsxfun(@times,nx,radii);
ry=bsxfun(@times,ny,radii);
rz=bsxfun(@times,nz,radii)*cline_para.zRatio;
sx=bsxfun(@plus,xs,rx);
sy=bsxfun(@plus,ys,ry);
sz=bsxfun(@plus,zs,rz);


        if ismatrix(im)
            
fx=interp2(Fx,sx,sy,'*linear');
fy=interp2(Fy,sx,sy,'*linear');
fz=interp2(Fz,sx,sy,'*linear');

        else
        
fx=interp3(Fx,sx,sy,sz,'*linear');
fy=interp3(Fy,sx,sy,sz,'*linear');
fz=interp3(Fz,sx,sy,sz,'*linear');
        end
fr=nx.*fx+ny.*fy+nz.*fz;
fr=nansum(fr,2);
fx=nansum(fx,2);
fy=nansum(fy,2);
fz=nansum(fz,2);
fs=[fx fy fz fr];
fs(isnan(fs)) = 0;

            
            %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
            % viscos force between frames on all nodes:        
xyzstart=xyzs(1,1:3);
xyzend=xyzs(end,1:3);
            % stretch or shrink the snake at the ends            
%             if cline_para.stretch_ends_flag
%                 s_head = xyzs(1,:) - xyzs(4, :);
%                 s_head = s_head/norm(s_head);
%                 fs(1, :) = fs(1, :) + s_head .* cline_para.stretching_force_factor(1);
%                 s_tail = xyzs(end, :) - xyzs(end - 3, :);
%                 s_tail = s_tail/norm(s_tail);
%                 fs(end, :) = fs(end, :) + s_tail .* cline_para.stretching_force_factor(2);
%             end
            %calculate the new position of snake
            xyzs(:,1:3) = Ainv * (cline_para.gamma*xyzs(:,1:3)  + cline_para.kappa*fs(:,1:3) );
            xyzs(:,4) = Arinv * (cline_para.gamma*xyzs(:,4)  + cline_para.kappaR*fs(:,4) );
            
            % take average value at the end of the cline_para.iterations
            if i > cline_para.iterations - 50
                xyzs_avg = xyzs_avg + xyzs;
                counter = counter + 1;
            end
%             xyzs(1,1:3)=xyzstart;
%             xyzs(end,1:3)=xyzend;
            L=length(xyzs)-1;
s=[0,cumsum(sqrt(diff(xyzs(:,1)).^2+diff(xyzs(:,2)).^2+diff(xyzs(:,3)).^2))']';
 
    xyzs(:,1)=interp1(s,xyzs(:,1),0:max(s)/L:max(s),'spline');
    xyzs(:,2)=interp1(s,xyzs(:,2),0:max(s)/L:max(s),'spline');
    xyzs(:,3)=interp1(s,xyzs(:,3),0:max(s)/L:max(s),'spline');
    xyzs(:,4)=interp1(s,xyzs(:,4),0:max(s)/L:max(s),'spline');
            
    if cline_para.show_flag==2
imagesc(Fx(:,:,26));
hold on
%scatter(xyzs(:,1),xyzs(:,2));
scatter(sx(:),sy(:));
hold off
            axis equal
   %         xlim([ 0 600]);ylim([ 0 600]);

drawnow
    end

                
%         figure(1);
% %        subplot(1, 2, 1);
%         cla;
%         imshow(stack_z_size_max(:,:,k),[]); 
%         hold on;
%         plot(xyzs(:,1), xyzs(:,2), 'r.-');
%         hold off;
%         disp(i);
%             
         end
         xyzs_avg = xyzs_avg / counter;
%         xyzs_all(:,:,k) = xyzs_avg;
% 
%         if k == 1 && pass == 1
%             xyzs_predicted(:,:,1) = xyzs_avg;
%             xyzs_predicted(:,:,2) = xyzs_avg;        
%         elseif k < length(time_idx) && pass == 1
%             xyzs_predicted(:,:,k+1) = xyzs_avg; 
%         end    

%         if pass == 2
%             figure(2)
%             cla;
%          imshow(stack_z_size_max(:,:,k),[]); 
%          hold on;
%          plot(xyzs(:,1), xyzs(:,2), 'r.-');
%          hold off;
%             writeVideo(vidObj,getframe);
%         end
%         figure(1); 
%         subplot(2, 1, 1);
%         cla;
%         imshow(stack_z_size_max(:,:,k),[]); 
%         hold on;
%         plot(xyzs_predicted(:,1, k), xyzs_predicted(:,2, k), 'c.-');
%         plot(xyzs(:,1), xyzs(:,2), 'r.-');
%         hold off;
%         subplot(2, 1, 2);
%         plot(xyzs(:, 1), xyzs(:, 3), '.')
%         axis([0, col, 0, 10]);
    
    % smooth predicted position for the 2nd pass
%     xyzs_all_pad = padarray(xyzs_all, [0 0 2], 'replicate');
%     xyzs_predicted = convn(xyzs_all_pad, reshape([1/5, 1/5, 1/5, 1/5, 1/5], 1,1,5), 'same');
%     xyzs_predicted = xyzs_predicted(:,:,3:end-2);
end
% close(vidObj);

% save data:
% save('xyzs_all.mat', 'xyzs_all');
% clearvars -except xyzs_all hyper_stack stack_z_size_max stack_t_size;