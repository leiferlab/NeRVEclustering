function [xyzs_avg,Is,EOut] = ActiveContourFit_wormRef4(im,tip_image, cline_para, cline_initial,refIdx,refPt)
% Open active contour program based on Yi's cell centerline fitting, which
% is an open version of the Kass 1988 active contour paper.
%Added several things to make the centerline work better for crappy worm
% images.

%%% INPUTS
%im  - input image, currently only 1D images are accepted.
%cline_para  - fitting parameters for active contour
%cline_initial - initial snake that is to be relaxed
%refIdx - index of point in snake that should be attached via spring to RefPt
%refPt - point which the refIdx should be attached to, based on low mag
%fluor image
% initialize snake
%    cline_initial=distanceInterp(cline_initial(2:end-1,:),100);

xyzs = distanceInterp(cline_initial(2:end-1,:),100);
[m, n ] = size(xyzs);

%intializing this value for later
ftail=0;
refL=cline_para.refL;
% for manually clicked tips
if isfield(cline_para,'head_pt')
    head_pt=cline_para.head_pt;
    tail_pt=cline_para.tail_pt;
    tip_flag=1;
else
    head_pt=[];
    tail_pt=[];
    tip_flag=0;
end

dIgnore=true(m,m);
dIgnore=triu(dIgnore,-20) & tril(dIgnore,20);
%dIgnore=[dIgnore rot90(dIgnore,2)];
repulsionD=cline_para.repulsionD;
repulsionD2=repulsionD*1.5;
gKernal2=gausswin(25,2);
gKernal2=gKernal2/sum(gKernal2);

[row ,col, stack_z_size] = size(im);
minDist=Inf;
minCounter=0;

% populating the penta diagonal matrix for active contour
if isfield(cline_para,'Ainv')
    Ainv=cline_para.Ainv;
else
    A = zeros(m,m);
    brow = zeros(1,m);
    brow(1,1:5) = [cline_para.CLbeta,...
        -(cline_para.CLalpha + 4*cline_para.CLbeta),...
        (2*cline_para.CLalpha + 6 *cline_para.CLbeta),...
        -(cline_para.CLalpha + 4*cline_para.CLbeta),...
        cline_para.CLbeta];
    
    A(1, 1:3) = [cline_para.CLbeta, -2*cline_para.CLbeta, cline_para.CLbeta];
    A(2, 1:4) = [-cline_para.CLalpha-2*cline_para.CLbeta,...
        2*cline_para.CLalpha+5*cline_para.CLbeta,...
        -cline_para.CLalpha-4*cline_para.CLbeta,...
        cline_para.CLbeta];
    for i=3:m-2
        A(i,:) = brow;
        brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
    end
    A(m-1, m-3:m) = [cline_para.CLbeta,...
        -cline_para.CLalpha-4*cline_para.CLbeta,...
        2*cline_para.CLalpha+5*cline_para.CLbeta,...
        -cline_para.CLalpha-2*cline_para.CLbeta];
    
    A(m, m-2:m) = [cline_para.CLbeta,...
        -2*cline_para.CLbeta,...
        cline_para.CLbeta];
    
    [L, U] = lu(A + cline_para.gamma .* eye(m,m));
    Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization
end

%% filter image to "normalize" intensity
% make smoothed image
gauss_kernel=fspecial('gaussian',25,5);
imSmooth=imfilter(im,gauss_kernel);

%testing better tip tracking
initial_head_I=1;%tip_image(round(refPt(2)),round(refPt(1)))/10;

%gradient
[fx, fy] = gradient(imSmooth*cline_para.gradient_force); %computing the gradient


%% prepare for cline_para.iterations
xyzs_avg = zeros(size(xyzs));
counter = 0;
%moving the snake in each iteration
for i=1:cline_para.iterations
    
    %% calculate image gradient force on the nodes
    
    s_ind = round(xyzs);
    f_image = zeros(m, n);
    in_image_flag = all((s_ind>0)&(bsxfun(@le, s_ind, [col, row])), 2);
    gradient_ind = sub2ind_nocheck([row, col], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1));
    f_image(in_image_flag, :) = [fx(gradient_ind), fy(gradient_ind), ];
    Is=zeros(size(f_image(:,1)));
    Is(in_image_flag)=tip_image(gradient_ind);
    if i==1
        Is_initial=Is;
        % for tail intenisty, take half the last 10 points, bounded by

            tail_I_init=mean(Is_initial(90:end))/2;
            tail_I_init=min(tail_I_init,.2);
            tail_I_init=max(tail_I_init,.06);
    end
    
    f_image=sqrt(abs(f_image)).*sign(f_image)*mean(abs(f_image(:)));
    fs=f_image;
 
    %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
    % viscos force between frames on all nodes:
    
    %% stretch or shrink the snake at the ends
    s=sqrt(sum(diff(xyzs).^2,2));
    
    xDiff=Gradient1(xyzs(:,1));
    yDiff=Gradient1(xyzs(:,2));
    xyDiff=[xDiff,yDiff];
    tVector=bsxfun(@rdivide,xyDiff,sqrt(sum(xyDiff.^2,2)));
    nVector_x=Gradient1(tVector(:,1));
    nVector_y=Gradient1(tVector(:,2));
    nVector=[nVector_x nVector_y];
    deltaS_all=s-refL;
    deltaS=mean(deltaS_all(refIdx:end));
    %     s_spring=-cline_para.springForce*...
    %         (f_forward+f_back);
    %
    
    %% pull on refernce point
    % pull on the reference point, pulling it towards the center found by
    % the lowmag fluor image. 
    refDistance=refPt-xyzs(refIdx,:);
    f_refSpring=refDistance*cline_para.refSpring;
    
    
    % bring down the ref spring force to tangent and normal componentes 
    f_refSpring_t=(f_refSpring*tVector(refIdx,:)')*tVector(refIdx,:);
    f_refSpring_n=f_refSpring-f_refSpring_t;
    
    % constrain the magnitude of the force from the spring by bounding the
    % distance, but only for the normal component
    
    abs_refDistance=sqrt(sum(refDistance.^2));
    abs_refDistance_n=abs_refDistance-10;
    abs_refDistance_n(abs_refDistance_n>60)=60;
    abs_refDistance_n(abs_refDistance_n<0)=0;
    f_refSpring_n=f_refSpring_n/abs_refDistance*abs_refDistance_n;
    f_refSpring_t=f_refSpring_t/abs_refDistance*(abs_refDistance);
    f_refSpring=f_refSpring_t+f_refSpring_n;

    %pull on all points from head to 2*refpt
    fs(1:2*refIdx,:)=bsxfun(@plus, fs(1:2*refIdx,:),f_refSpring);
    
    s_spring=nVector*deltaS*cline_para.refSpring*1;
    if i>cline_para.iterations/2
        s_spring=conv2(s_spring,gKernal2,'same');
    end
    fs(2:end-1,:)=fs(2:end-1,:)+s_spring(2:end-1,:);
    
    
    %% tip repulsion
    
    %The body of the snake repels the tips if they are close by. 
    % Characteristic distances is repulsionD
    
    %make distance matrices, using only points near the tip to rest of the
    %body, x and y seperately then together
    xmat=bsxfun(@minus,xyzs(:,1)',xyzs(:,1));
    ymat=bsxfun(@minus,xyzs(:,2)',xyzs(:,2));
    dmat=sqrt(xmat.^2+ymat.^2);
    
    %say ignored points are far
    dmat(dIgnore)=100;
    
    %get vectors between points to project the forces along, col i is all
    %vectors pointing towards point i,
    v_x=xmat./dmat;
    v_y=ymat./dmat;
    f_end=zeros(size(dmat));
    %only apply forces to the close points
    close_pts=dmat<repulsionD*2.5;
    % for of repulsion potential, modelled by soft shell
    if any(close_pts(:))
        f_end(close_pts)=repulsionD^6./(dmat(close_pts).^6+repulsionD^6);
    end
    
    %sum the forces to points and the end
    f_repulsion_x=sum(f_end.*v_x);
    f_repulsion_y=sum(f_end.*v_y);
    f_repulsion=cline_para.endRepulsion*...
        [f_repulsion_x' f_repulsion_y'];

    fs=fs+f_repulsion;
    
    %% repel tail and body from refPt
    refPt_v=bsxfun(@minus,xyzs(50:end,:),refPt);
    refPt_d=sqrt(sum(refPt_v.^2,2));
    f_end_ref=zeros(length(refPt_d),2);
    close_pts=refPt_d<repulsionD2*2.5;
    if any(close_pts(:))
        ref_repulsion=repulsionD2^6./(refPt_d.^6+repulsionD2^6);
        f_end_ref=bsxfun(@times , ref_repulsion./refPt_d,refPt_v);
        
    end
    f_end_ref=cline_para.endRepulsion*f_end_ref*5;
    fs(50:end,:)=fs(50:end,:)+f_end_ref;
    
    %% stretch ends based on intensity, still buggy
    if cline_para.stretch_ends_flag
        if i>5
            %if intensity of head point is above zero, grow if the current
            %intensity is less than the inital head point intensity, add
            %growing force
            if  Is(1)>0
                fhead=cline_para.endkappa*(Is(1)-initial_head_I-4*deltaS.^2*sign(deltaS));
                fhead(fhead>15)=15;
            else
                %if initial head intensity is zero, include stretch force 
                %based on length (deltaS). 
                fhead=nnz(Is(1:10)==0)*cline_para.endkappa*(-initial_head_I)-deltaS;
            end
            %fo tail stretching force, check the tail intensity using a
            %number of different end ranges. Seems to be working so far. 
            tail_I_20=mean(Is(end-20:end));
            tail_I_10=mean(Is(end-10:end));
            tail_I_5=mean(Is(end-5:end));
            

            ftail=cline_para.endkappa*(tail_I_5-1-deltaS.^2.*sign(deltaS));
            
            ftail(ftail>5)=5;
            ftail(ftail<-5)=-5;
        else
            ftail=0;
            fhead=0;
        end
        % stretch the head along s_head (head direction)
        s_head = xyzs(1,:) - xyzs(5, :);
        s_head = s_head/norm(s_head);
        fs(1, :) =fs(1, :)+  s_head .* fhead*cline_para.stretching_force_factor(1);
        s_tail = xyzs(end, :) - xyzs(end - 4, :);
        
        % stretch the tail along s_tail
        s_tail = s_tail/norm(s_tail);
        fs(end, :) =  fs(end, :)+ s_tail .* ftail*cline_para.stretching_force_factor(2);
        % add random noise, temperature decreases over time. 
        fhead_heat=cline_para.heat*(rand(1,2)-.5).*exp(-i/cline_para.iterations*10);
        ftail_heat=cline_para.heat*(rand(1,2)-.5).*exp(-i/cline_para.iterations*10);
        fs(1,:)=fs(1,:)+fhead_heat;
        fs(end,:)=fs(end,:)+ftail_heat;
    end
    
    %% memory force pushes the centerline towards the initial centerline.
    f_memory=-xyzs+cline_initial;
    %scale that force by the worm intensity. 
    if any(Is_initial) %normalization breaks if all Is_initial is zero
        f_memory=bsxfun(@times, f_memory, Is_initial./mean(Is_initial));
        f_memory=f_memory*cline_para.memForce;
        fs=fs+f_memory;
    end
    %% plot some results if show is on
    
    if cline_para.showFlag && ~mod(i,cline_para.showFlag)
        imagesc(imSmooth(:,:,1));colorbar
        hold on
        plot(xyzs(:,1),xyzs(:,2),'r')
        scatter(refPt(1),refPt(2),'gx');
        plot([xyzs(refIdx,1) refPt(1)],[xyzs(refIdx,2) refPt(2)],'g');
        plotfs=fs;
        plotfs(refIdx,:)=0;
        quiver(xyzs(:,1),xyzs(:,2),plotfs(:,1),plotfs(:,2),'black');
        plot(xyzs(:,1),xyzs(:,2))
        hold off
        drawnow
    end
    %% move points
    
    %calculate the new position of snake
    xyzsNew = Ainv * (cline_para.gamma*xyzs + cline_para.kappa*fs);
    % take average value at the end of the cline_para.iterations
    deltaxyzs=xyzsNew-xyzs;
    xyzs=xyzsNew;
    
    %enforce tip locations if manually clicked
    if any(head_pt)
        xyzs(1,:)=head_pt;
        xyzs(end,:)=tail_pt;
    end
    
    distanceTravelled=sum(sqrt(sum(diff(deltaxyzs).^2,2)));
    if i > cline_para.iterations - 50
        %turn of ref spring to try to get the contour to relax smoothly
        %  cline_para.refSpring=0;
        xyzs_avg = xyzs_avg + xyzs;
        counter = counter + 1;
    end
    if distanceTravelled<minDist || i<cline_para.iterations/2
        minCounter=0;
        minDist=distanceTravelled;
    else
        minCounter=minCounter+1;
    end
    
    %% every 10 frames, remove 4 points from head and tail and try to have CL regrow
    %trying to stop tips from getting lost
    if i > cline_para.iterations - 50 && ~mod(i,20) &&  cline_para.stretch_ends_flag
        xyzs=distanceInterp(xyzs(4:end-3,:),100);
    end
    
    %% Try to fix places where worm crosses itself and try to correct it
    cross_bool= doesCross(xyzs);
    if cross_bool
        display('Cross detected')
        %find linesegments that cross
        % cross_1(i) and cross_2(i) are the indeces of the ith cross. 
        [~,cross_1,cross_2]= doesCross(xyzs);
        cross_list=[];
        for iCross=1:length(cross_1)
            subCrossList=cross_1(iCross):cross_2(iCross);
            %if its a small loop, excise the loop. 
            if length(subCrossList)<25
                cross_list=[cross_list subCrossList];
                %if its a small loop at the head, cut off part of the head
            elseif cross_1(iCross)<3
                cross_list=[cross_list 1:3];
            end
        end
        xyzs(cross_list,:)=[];
        xyzs=distanceInterp(xyzs,100);
    elseif ~mod(i,20) 
        xyzs=distanceInterp(xyzs,100);
    end
    
end

%do some averageing at the end to smooth out the snake.
xyzs_avg = xyzs_avg / counter;
xyzs_avg=distanceInterp(xyzs_avg,100);
if ftail<-.1
    xyzs_avg=distanceInterp(xyzs_avg(1:97,:),100);
end
%caluclate the internal energy of the snake for output
EOut=1/2*(cline_para.CLalpha*(sum(sum(diff(xyzs_avg,1,1).^2))) + ...
    cline_para.CLbeta*(sum(sum(diff(xyzs_avg,2,1).^2))));


end

function g=Gradient1(f)
%faster 1D gradient function
n=length(f);
h=(1:n)';
g(1,:) = (f(2) - f(1))/(h(2)-h(1));
g(n,:) = (f(n) - f(n-1))/(h(end)-h(end-1));

h2 = h(3:n) - h(1:n-2);
g(2:n-1) =(f(3:n)-f(1:n-2))./h2;
end
