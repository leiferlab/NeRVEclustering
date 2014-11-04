function xyzs_avg = ActiveContourFit2_1D_2(im, cline_para, cline_initial)

if isempty(cline_para);
    cline_para.iterations=100;
cline_para.threshold=.1;
cline_para.stiff=0.05; %.0005
cline_para.stiff3d=1; %.1
cline_para.alpha =10 ;
cline_para.beta = 80000;

cline_para.inter_frame_viscos_factor=0; %0
cline_para.kappa=500;
cline_para.gamma=1; %step size
cline_para.gradflag=0;
cline_para.gradient_force =.5;%5
cline_para.interpolate_gradient_flag =0;
cline_para.show_flag=1; %1 to show end of every frame, 2 to show every compelted slice
cline_para.movie_flag=0;
end

if isempty(cline_initial)
    [I,centerLine]=nanmax((im'));
    I(isnan(I))=0;
%    centerLine(centerLine<(size(im,2)*.1)| centerLine>(size(im,2)*.9))=size(im,2)/2;
%    centerLine=conv(centerLine.*I,ones(1,300),'same')./conv(I,ones(1,300),'same');
    centerLine(isnan(centerLine))=size(im,2)/2;
    
    cline_initial=[centerLine',(1:size(im,1))'];
end





[row col, stack_z_size] = size(im);

% initialize snake
xyzs = cline_initial;
%s_interp = 0:0.1:s(end);
%xyzs = interp1(s, cline_initial, s_interp, 'spline');


[m n] = size(xyzs);
    
% populating the penta diagonal matrix
A = zeros(m,m);
brow = zeros(1,m);
brow(1,1:5) = [cline_para.beta, -(cline_para.alpha + 4*cline_para.beta), (2*cline_para.alpha + 6 *cline_para.beta), -(cline_para.alpha + 4*cline_para.beta), cline_para.beta];
A(1, 1:3) = [cline_para.beta, -2*cline_para.beta, cline_para.beta];
A(2, 1:4) = [-cline_para.alpha-2*cline_para.beta, 2*cline_para.alpha+5*cline_para.beta, -cline_para.alpha-4*cline_para.beta, cline_para.beta];
for i=3:m-2
    A(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end
A(m-1, m-3:m) = [cline_para.beta, -cline_para.alpha-4*cline_para.beta, 2*cline_para.alpha+5*cline_para.beta, -cline_para.alpha-2*cline_para.beta];
A(m, m-2:m) = [cline_para.beta, -2*cline_para.beta, cline_para.beta];
[L U] = lu(A + cline_para.gamma .* eye(m,m));
Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization




        %% gradient
        
        [fy, fx] = gradient(im*cline_para.gradient_force); %computing the gradient
        fx=zeros(size(fx));
%         fy=fy;
if cline_para.show_flag
         figure;
 imagesc(im);
hold on
end
        %% prepare for cline_para.iterations
        xyzs_avg = zeros(size(xyzs));
        counter = 0;
        %moving the snake in each iteration
        for i=1:cline_para.iterations;
            %calculate image gradient force on the nodes
             % actual interpolation, slow!!!
                xs=xyzs(:,1);
                ys=xyzs(:,2);
                fs = [interp2(fy,xs,ys,'*linear'),...
                    interp2(fx,xs,ys,'*linear')]     ;     
                
                
                fs(isnan(fs)) = 0;
%                 if i<cline_para.iterations/2
%                 fs(end,1)=(xs(round(length(xs)*2/3))-xs(end))*.01;
%                 fs(1,1)=(xs(round(length(xs)*1/3))-xs(1))*.01;
%                 end

%  fs(end,:)=5*nanmean(fs(round(length(fs)/3):end,:));
%   fs(1,:)=1*nanmean(fs(1:round(length(fs)/3),:));

            %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
            % viscos force between frames on all nodes:        
xyzstart=xyzs(1,:);
xyzend=xyzs(end,:);
            % stretch or shrink txhe snake at the ends            
%             if cline_para.stretch_ends_flag
%                 s_head = xyzs(1,:) - xyzs(4, :);
%                 s_head = s_head/norm(s_head);
%                 fs(1, :) = fs(1, :) + s_head .* cline_para.stretching_force_factor(1);
%                 s_tail = xyzs(end, :) - xyzs(end - 3, :);
%                 s_tail = s_tail/norm(s_tail);
%                 fs(end, :) = fs(end, :) + s_tail .* cline_para.stretching_force_factor(2);
%             end
            %calculate the new position of snake
            xyzs(:,1) = Ainv * (cline_para.gamma*xyzs(:,1) + cline_para.kappa*fs(:,1));
            % take average value at the end of the cline_para.iterations
            if i > cline_para.iterations - 20
                xyzs_avg = xyzs_avg + xyzs;
                counter = counter + 1;
            end
            
            %end points are fixed

            xyzs(1,2)=xyzstart(2); 
            xyzs(end,2)=xyzend(2);
            L=length(xyzs)-1;
            %redistribute points evenly
s=[0,cumsum(sqrt(diff(xyzs(:,1)).^2+diff(xyzs(:,2)).^2))']';
 
    xyzs(:,1)=interp1(s,xyzs(:,1),0:max(s)/L:max(s),'spline');
    xyzs(:,2)=interp1(s,xyzs(:,2),0:max(s)/L:max(s),'spline');
          
    %view contour
    if cline_para.show_flag
    if exist('h','var')
    set(h,'Xdata',xyzs(:,1),'Ydata',xyzs(:,2));
    else
    h=scatter(xyzs(:,1),xyzs(:,2));
    axis equal
    end
    drawnow
    end
%     
          
        end
         %avg last 50 iterations
         xyzs_avg = xyzs_avg / counter;

end
