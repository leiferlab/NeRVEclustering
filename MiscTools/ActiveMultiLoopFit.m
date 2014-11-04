function countourOut = ActiveMultiLoopFit(hyper_stack, cline_para, cline_initial)
imDims=2;%ndims(hyper_stack);

A=[];
% cline_para.alpha =10; 
% cline_para.beta = 200;
nLoops=length(cline_initial);
[row col, stack_z_size] = size(hyper_stack);
for iLoops=1:nLoops

% initialize snake
xyzs_i = cline_initial(iLoops).xyzs;

%s_interp = 0:0.1:s(end);
%xyzs = interp1(s, cline_initial, s_interp, 'spline');
[m n] = size(xyzs_i);
if iLoops==1
    startIdx(iLoops)=1;
    endIdx(iLoops)=m;
else
    startIdx(iLoops)=endIdx(iLoops-1)+1;
    endIdx(iLoops)=endIdx(iLoops-1)+m;
end

    
   
% populating the penta diagonal matrix
subA = zeros(m,m);
brow = zeros(1,m);
brow(1,1:5) = [cline_para.CLbeta, -(cline_para.CLalpha + 4*cline_para.CLbeta), (2*cline_para.CLalpha + 6 *cline_para.CLbeta), -(cline_para.CLalpha + 4*cline_para.CLbeta), cline_para.CLbeta];
brow=circshift(brow,[0,-2]);

for i=1:m
    subA(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end

A=blkdiag(A,subA);
end
xyzs=cell2mat({cline_initial.xyzs}');
if size(xyzs,2)==2;
    xyzs=[xyzs,ones(length(xyzs),1)];
end


[L U] = lu(A + cline_para.gamma .* eye(length(A)));

Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization

% find the center line of the cell

im = zeros(row, col, stack_z_size);
% two pass algorithm

% vidObj = VideoWriter('tmp.avi', 'Motion Jpeg AVI');
% vidObj.FrameRate = 15;
% vidObj.Quality = 95;
% open(vidObj);
    %k: stack_z_size index
        %obtain stack_z_size(t)
for kSlice=1:size(hyper_stack,3);
        %% filter image to "normalize" intensity
im=hyper_stack(:,:,kSlice);


if cline_para.gradFlag
                kernal=fspecial('sobel');
            kernal3(:,:,1)=kernal;

            kernal3r=kernal;
            kernal3c=kernal';
            hostack=convn(im,kernal3r,'same');
            vertstack=convn(im,kernal3c,'same');
            
            im=sqrt(hostack.^2+vertstack.^2);
end

        %gradient
        if ndims(im)==3;
        [fx, fy, fz] = gradient(im); %computing the gradient
        else
         [fx, fy] = gradient(im); %computing the gradient
fz=zeros(size(fx));
        end
        %% prepare for cline_para.iterations
        xyzs_avg = zeros(size(xyzs));
        counter = 0;
        %moving the snake in each iteration
        for i=1:cline_para.iterations;
            %calculate image gradient force on the nodes
            if 1                % actual interpolation, slow!!!
                xs=xyzs(:,1);
                ys=xyzs(:,2);
                zs=xyzs(:,3);
                if imDims==3;
                fs = [interp3(fx,ys,xs, zs), interp3(fy,ys,xs, zs), interp3(fz,ys,xs, zs)];
                else 
                fs = [interp2(fx,xs,ys), interp2(fy,xs,ys), interp2(fz,ys,xs)];
                end
                fs(isnan(fs)) = 0;
            else
                s_ind = round(xyzs);        
                fs = zeros(m, 3);    
                in_image_flag = all((s_ind>0)&(bsxfun(@le, s_ind, [col, row, stack_z_size])), 2);
                gradient_ind = sub2ind([row, col, stack_z_size], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1), s_ind(in_image_flag, 3));
                fs(in_image_flag, :) = -[fx(gradient_ind), fy(gradient_ind), fz(gradient_ind)];
            end
            
            %% repulsion force
        Frepulsion=zeros(size(xyzs));
        for ipoly=2:nLoops
            for jpoly=1:ipoly-1
                xmat=bsxfun(@minus,xyzs(startIdx(ipoly):endIdx(ipoly),1)',...
                    xyzs(startIdx(jpoly):endIdx(jpoly),1));
                ymat=bsxfun(@minus,xyzs(startIdx(ipoly):endIdx(ipoly),2)',...
                    xyzs(startIdx(jpoly):endIdx(jpoly),2));
                zmat=bsxfun(@minus,xyzs(startIdx(ipoly):endIdx(ipoly),2)',...
                    xyzs(startIdx(jpoly):endIdx(jpoly),2));
                
                dmat=sqrt(xmat.^2+ymat.^2+zmat.^2);
                xmat=xmat./dmat;
                ymat=ymat./dmat;
                zmat=zmat./dmat;
                
                Fmat=(1-dmat/cline_para.objsize).^2;
                Fmat(dmat>cline_para.objsize)=0;
                
                Frepulsionx=Fmat.*xmat;
                Frepulsiony=Fmat.*ymat;
                Frepulsionz=Fmat.*zmat;
                Frepulsion(startIdx(ipoly):endIdx(ipoly),:)=...
                    Frepulsion(startIdx(ipoly):endIdx(ipoly),:)-...
                    [sum(Frepulsionx,1)',sum(Frepulsiony,1)',sum(Frepulsionz,1)'];
                
                Frepulsion(startIdx(jpoly):endIdx(jpoly),:)=...
                    Frepulsion(startIdx(jpoly):endIdx(jpoly),:)+...
                    [sum(Frepulsionx,2),sum(Frepulsiony,2),sum(Frepulsionz,2)];
                
                
                
            end
        end
        fs=fs+Frepulsion*cline_para.repulsion;
            
        %% loop spring force
 
        
        loopCM=zeros(length(startIdx),3);
        for iLoops=1:nLoops
            loopCM(iLoops,:)=mean(xyzs(startIdx(iLoops):endIdx(iLoops),:),1);
        end
        distanceVector=normr(diff(loopCM,[],1));
        
        
        loopDistance=sqrt(sum(diff(loopCM,[],1).^2,2));
        loopDelta=loopDistance-cline_para.distance0';
        loopForce=loopDelta*cline_para.loopSpringForce;
        
        loopForce=bsxfun(@times,loopForce,distanceVector);
        for iLoops=1:nLoops
            if iLoops<nLoops
            fs(startIdx(iLoops):endIdx(iLoops),:)=...
                bsxfun(@plus,fs(startIdx(iLoops):endIdx(iLoops),:),loopForce(iLoops,:));
            end
            if iLoops>1
            fs(startIdx(iLoops-1):endIdx(iLoops-1),:)=...
                bsxfun(@minus,fs(startIdx(iLoops-1):endIdx(iLoops-1),:),loopForce(iLoops-1,:));
            end
        end
        
                    
        
        
            %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
            % viscos force between frames on all nodes:        

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
            xyzsnew = Ainv * (cline_para.gamma*xyzs + cline_para.kappa*fs);
            % take average value at the end of the cline_para.iterations
           deltaxyzs=xyzsnew-xyzs;
           if i<cline_para.groupIterations;
               deltaxyzs=mean(deltaxyzs);
               xyzs=bsxfun(@plus, xyzs,deltaxyzs);
           else
               xyzs=xyzs+deltaxyzs;
           end
           
            
            if i > cline_para.iterations - 50
                xyzs_avg = xyzs_avg + xyzs;
                counter = counter + 1;
            end

            for iLoops=1:nLoops
                subLoop=xyzs(startIdx(iLoops):endIdx(iLoops),:);
                subLoop=[subLoop;subLoop(1,:)];
            L=length(subLoop)-1;
s=[0,cumsum(sqrt(diff(subLoop(:,1)).^2 ...
+diff(subLoop(:,2)).^2 ...
+diff(subLoop(:,3)).^2))']';
 
    subLoop(:,1)=interp1(s,subLoop(:,1), ...
        0:max(s)/L:max(s),'spline');
   subLoop(:,2)=interp1(s,subLoop(:,2), ...
        0:max(s)/L:max(s),'spline');
    subLoop(:,3)=interp1(s,subLoop(:,3), ...
        0:max(s)/L:max(s),'spline');
    
    subLoop=subLoop(1:end-1,:);
    xyzs(startIdx(iLoops):endIdx(iLoops),:)=subLoop;
            end
    if cline_para.show_flag>1
        if mod(i,cline_para.show_flag)==0
        if imDims==3
             scatter3(xyzs(:,1),xyzs(:,2),xyzs(:,3));
            axis equal
pause(.01)
        else
         if ~exist('hfigure','var');
            hfigure=imagesc(hyper_stack(:,:,kSlice));axis equal;
            hold on
         else
             if ishandle(hfigure);delete(hfigure);end
             hfigure=imagesc(hyper_stack(:,:,kSlice));axis equal;
             hold on
         end
            if exist('hfigure2','var')
                if ishandle(hfigure2)
                    delete(hfigure2)
                    delete(hfigure3)
                end
            end
            title(num2str(kSlice))
            hfigure2=plot(xyzs(:,1),xyzs(:,2),'r');
            %scatter(xyzs(:,1),xyzs(:,2),'r');
            hfigure3=quiver(xyzs(:,1),xyzs(:,2),fs(:,1),fs(:,2));
            imEnergy=sum(interp2(im,xs,ys,'*cubic'));

            pause(.001)
        end
        
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
         
         for iLoops=1:nLoops
                subLoop=xyzs(startIdx(iLoops):endIdx(iLoops),:);
  
         countourOut(kSlice).cline(iLoops).xyzs=subLoop;
         end
         
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
        if cline_para.show_flag==1;
         if ~exist('hfigure','var');
            hfigure=imagesc(hyper_stack(:,:,kSlice));axis equal;
            hold on
         else
             if ishandle(hfigure);delete(hfigure);end
             hfigure=imagesc(hyper_stack(:,:,kSlice));axis equal;
         end
         
            if exist('hfigure2','var')
                if ishandle(hfigure2)
                    delete(hfigure2)
                end
            end
            hfigure2=plot(xyzs(:,1),xyzs(:,2),'r');
            %scatter(xyzs(:,1),xyzs(:,2),'r');
            imEnergy=sum(interp2(im,xs,ys,'*cubic'));
            
title(num2str(kSlice))
drawnow
        end
end
% close(vidObj);

% save data:
% save('xyzs_all.mat', 'xyzs_all');
% clearvars -except xyzs_all hyper_stack stack_z_size_max stack_t_size;