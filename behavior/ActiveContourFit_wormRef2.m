function [xyzs_avg,Is,EOut] = ActiveContourFit_wormRef2(hyper_stack, cline_para, cline_initial,refIdx,refPt)
% Open active contour program based on Yi's cell centerline fitting. Added
% several things to make the centerline work better for crappy worm
% images.

%%% INPUTS
%hyper_stack  - input image, currently only 1D images are accepted. 
%cline_para  - fitting parameters for active contour
%cline_initial - initial snake that is to be relaxed
%refIdx - index of point in snake that should be attached via spring to RefPt
%refPt - point which the refIdx should be attached to, based on low mag
%fluor image

refL=cline_para.refL;
gaussFilter=fspecial('gaussian',30,64);
    endR=cline_para.tipRegion;
    repulsionD=cline_para.repulsionD;


[row col, stack_z_size] = size(hyper_stack);
minDist=Inf;
minCounter=0;
% initialize snake
xyzs = cline_initial;
[m, n ] = size(xyzs);

% populating the penta diagonal matrix
if isfield(cline_para,'Ainv')
    Ainv=cline_para.Ainv;
else
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
end
% find the center line of the cell

im = zeros(row, col, stack_z_size);
% two pass algorithm


%% filter image to "normalize" intensity
im=hyper_stack;
imSmooth=filter2(gaussFilter,im,'same');

%gradient
if n==3
    [fy, fx, fz] = gradient(im*cline_para.gradient_force); %computing the gradient
elseif n==2
    [fx, fy] = gradient(im*cline_para.gradient_force); %computing the gradient
end

%% prepare for cline_para.iterations
xyzs_avg = zeros(size(xyzs));
counter = 0;
%moving the snake in each iteration
for i=1:cline_para.iterations;
    %calculate image gradient force on the nodes
    if 0                % actual interpolation, slow!!!
        xs=xyzs(:,1);
        ys=xyzs(:,2);
        if n==3
            zs=xyzs(:,3);
            fs = [interp3(fx,ys,xs, zs), interp3(fy,ys,xs, zs), interp3(fz,ys,xs, zs)];
        elseif n==2
            fs = [interp2(fx,xs,ys), interp2(fy,xs,ys)];
        end
        Is=interp(imSmooth,xs,ys)
        fs(isnan(fs)) = 0;
    else
        s_ind = round(xyzs);
        fs = zeros(m, n);
        in_image_flag = all((s_ind>0)&(bsxfun(@le, s_ind, [col, row])), 2);
        gradient_ind = sub2ind([row, col], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1));
        fs(in_image_flag, :) = [fx(gradient_ind), fy(gradient_ind), ];
        Is=zeros(size(fs(:,1)));
        Is(in_image_flag)=imSmooth(gradient_ind);
    end
    %fs=bsxfun(@plus,fs,(refPt-xyzs(refIdx,:))*cline_para.refSpring);
    %xyzs=distanceInterp(xyzs,length(xyzs));
    %[~,newRefIdx]=pdist2(xyzs,xyzs(refIdx,:),'euclidean','Smallest',1);

    %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
    % viscos force between frames on all nodes:
    f_memory=-xyzs+cline_initial;
    fs=fs+f_memory*cline_para.memForce;
    % xyzstart=xyzs(1,:);
    % xyzend=xyzs(end,:);
    % stretch or shrink the snake at the ends
    s=sqrt(sum(diff(xyzs).^2,2));
    
    xyDiff=diff(xyzs(:,:));
    tVector=bsxfun(@rdivide,xyDiff,s);
    deltaS=s-refL;
    deltaS_head=mean(deltaS(1:refIdx));
    deltaS=mean(deltaS);
    f_spring=bsxfun(@times,tVector,deltaS_head);
    f_forward=[0 0;f_spring];
    f_back=[-f_spring; 0  0];
    s_spring=-cline_para.springForce*...
        (f_forward+f_back);
    
        newRefIdx=12;
        refDistance=refPt-xyzs(newRefIdx,:);
            f_refSpring=refDistance*cline_para.refSpring;


    f_refSpring_t=(f_refSpring*tVector(refIdx,:)')*tVector(refIdx,:);
    f_refSpring_n=f_refSpring-f_refSpring_t;
    abs_refDistance=sqrt(sum(refDistance.^2));
    f_refSpring_n=f_refSpring_n/abs_refDistance*median([20, abs_refDistance-10,0]);
    f_refSpring_t=f_refSpring_t/abs_refDistance*min(abs_refDistance,10);

    fs(refIdx,:)=fs(refIdx,:)+f_refSpring_t+f_refSpring_n;
    
    
    %      s_spring(1,:)=tVector(2,:)*deltaS;
    %      s_spring(end,:)=-tVector(end-1,:)*deltaS;
    %      s_spring(2:end-1,:)=0;
    %       s_spring=(xyzs(1,:)-xyzs(refIdx,:));
    %       s_spring=0.1*(normr(s_spring)*(refIdx)*medS-s_spring);
    s_spring=s_spring(1:refIdx-1,:);
    fs(1:size(s_spring,1),:)=fs(1:size(s_spring,1),:)+s_spring;
    
        %tip repulsion
     xmat=bsxfun(@minus,xyzs([1:endR end-endR+1:end],1)',xyzs(:,1));
    ymat=bsxfun(@minus,xyzs([1:endR end-endR+1:end],2)',xyzs(:,2));

    dmat=sqrt(xmat.^2+ymat.^2);
    dmat(dmat==0)=100;
    v_x=xmat./dmat;
    v_y=ymat./dmat;
 
    f_end=repulsionD^6./(dmat.^6+repulsionD^6);
    endR2=round(1.5*endR);
    f_end(end-endR2:end,endR+1:end)=0;
    f_end(1:endR2,1:endR)=0;
    
f_tip_end_x=sum(f_end(:,endR+1:end)'*v_x(:,endR+1:end),2);
f_tip_end_y=sum(f_end(:,endR+1:end)'*v_y(:,endR+1:end),2) ;
f_tip_end=cline_para.endRepulsion*[f_tip_end_x f_tip_end_y];
    
f_tip_start_x=sum(f_end(:,1:endR)'*v_x(:,1:endR),2);
f_tip_start_y=sum(f_end(:,1:endR)'*v_y(:,1:endR),2) ;
f_tip_start=cline_para.endRepulsion*[f_tip_start_x f_tip_start_y];

% make tip forces tangent to CL
% f_tip_start=bsxfun(@times,f_tip_start,dot(tVector(1:endR,:),f_tip_start,2));
% f_tip_end=bsxfun(@times,f_tip_end,dot(tVector(1:endR,:),f_tip_end,2));




    if cline_para.stretch_ends_flag
        fhead=cline_para.endkappa*(Is(1)-mean(Is)*.3)-deltaS_head*1;
        ftail=cline_para.endkappa*(Is(end)-mean(Is)*.3)-deltaS*1;
        if fhead<0;
            fhead=fhead*4;
        end
        if ftail<0;
            ftail=ftail*4;
        end
        
        
        s_head = xyzs(1,:) - xyzs(4, :);
        s_head = s_head/norm(s_head);
        fs(1, :) = fs(1, :) + s_head .* fhead*cline_para.stretching_force_factor(1);
        s_tail = xyzs(end, :) - xyzs(end - 3, :);
        s_tail = s_tail/norm(s_tail);
        fs(end, :) = fs(end, :) + s_tail .* ftail*cline_para.stretching_force_factor(2);
        fhead_heat=cline_para.heat*(rand(1,2)-.5).*exp(-i/cline_para.iterations*10);
        ftail_heat=cline_para.heat*(rand(1,2)-.5).*exp(-i/cline_para.iterations*4);
        fs(end-2*refIdx:end, :)=fs(end-2*refIdx:end,:)*(1-exp(-i/cline_para.iterations*10));
        fs(1:refIdx, :)=fs(1:refIdx,:)*(1-exp(-i/cline_para.iterations*10));

        fs(1,:)=fs(1,:)+fhead_heat;
        fs(end,:)=fs(end,:)+ftail_heat;
    end
    fs(1:endR,:)=fs(1:endR,:)+f_tip_start;
 fs(end-endR+1:end,:)=fs(end-endR+1:end,:)+f_tip_end;
   

    %calculate the new position of snake
    xyzsNew = Ainv * (cline_para.gamma*xyzs + cline_para.kappa*fs);
    % take average value at the end of the cline_para.iterations
    deltaxyzs=xyzsNew-xyzs;
    xyzs=xyzsNew;
    distanceTravelled=sum(sqrt(sum(diff(deltaxyzs).^2,2)));
    if i > cline_para.iterations - 50
        
        %turn of ref spring to try to get the contour to relax smoothly
        cline_para.refSpring=0;
        xyzs_avg = xyzs_avg + xyzs;
        counter = counter + 1;
    end
    %             xyzs(1,:)=xyzstart;
    %             xyzs(end,:)=xyzend;
    % if ~mod(i,15)
    % 
    % end
    if distanceTravelled<minDist || i<cline_para.iterations/2
        minCounter=0;
        minDist=distanceTravelled;
    else
        minCounter=minCounter+1;
                cline_para.refSpring=0;
                xyzs=distanceInterp(xyzs,length(xyzs));

    end
    
    if minCounter>50;
        xyzs_avg=xyzs;
        counter=1;
        break;
    end
    if doesCross(xyzs(1:30,:));
        xyzs=distanceInterp(xyzs(refIdx:end,:),length(xyzs));
    end
    
    
    if cline_para.showFlag
        if ~mod(i,cline_para.showFlag)
            imagesc(fx(:,:,1));
            hold on
            plot(xyzs(:,1),xyzs(:,2),'r')
            scatter(refPt(1),refPt(2),'gx');
            plot([xyzs(newRefIdx,1) refPt(1)],[xyzs(newRefIdx,2) refPt(2)],'g');
            quiver(xyzs(:,1),xyzs(:,2),fs(:,1),fs(:,2),'black');
            hold off
            drawnow
            
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
%     %
%     EOut=1/2*(cline_para.CLalpha*(sum(sum(diff(xyzs,1,2).^2))) + ...
%     cline_para.CLbeta*(sum(sum(diff(xyzs,2,2).^2))))
end
xyzs_avg = xyzs_avg / counter;
EOut=1/2*(cline_para.CLalpha*(sum(sum(diff(xyzs_avg,1,1).^2))) + ...
    cline_para.CLbeta*(sum(sum(diff(xyzs_avg,2,1).^2))));
% xyzs_avg=distanceInterp(xyzs_avg,length(xyzs_avg));
