function xyzs_avg = ActiveContourFit2_psudo2D(im, cline_para, cline_initial)

if isempty(cline_para);
    cline_para.iterations=100;
cline_para.threshold=.1;
cline_para.alpha =100 ;
cline_para.beta = 10000;
cline_para.alpha3D=.0001;
cline_para.beta3D=1000;

cline_para.inter_frame_viscos_factor=0; %0
cline_para.kappa=800;
cline_para.gamma=1; %step size
cline_para.gradflag=0;
cline_para.gradient_force =.5;%5
cline_para.interpolate_gradient_flag =0;
cline_para.show_flag=0; %1 to show end of every frame, 2 to show every compelted slice
cline_para.movie_flag=0;
end

if isempty(cline_initial)
    [I,centerLine]=nanmax(im,[],2);
    I(isnan(I))=0;
    I(I<quantile(I(:),.3))=0;
%    centerLine(centerLine<(size(im,2)*.1)| centerLine>(size(im,2)*.9))=size(im,2)/2;
%    centerLine=conv(centerLine.*I,ones(1,300),'same')./conv(I,ones(1,300),'same');
    centerLine(isnan(centerLine))=size(im,2)/2;
    [centerLineSpacex,centerLineSpacey]=meshgrid(1:size(centerLine,3),1:size(centerLine,1));
    sCenterLine=squeeze(centerLine);
    f=fit([centerLineSpacex(:),centerLineSpacey(:)],sCenterLine(:),'poly22', 'Weight',(I(:)));
    fCenterLine=f(centerLineSpacex,centerLineSpacey);
    centerLine=reshape(fCenterLine,size(centerLine));
    cline_initial=cat(2,centerLine,repmat((1:size(im,1))',1,1,size(centerLine,3)));


end

%%



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





%% make time pentadiagmatrix

B= zeros(stack_z_size,stack_z_size);
brow = zeros(1,stack_z_size);
brow(1,1:5) = [cline_para.beta3D, -(cline_para.alpha3D + 4*cline_para.beta3D), (2*cline_para.alpha3D + 6 *cline_para.beta3D), -(cline_para.alpha3D + 4*cline_para.beta3D), cline_para.beta3D];
B(1, 1:3) = [cline_para.beta3D, -2*cline_para.beta3D, cline_para.beta3D];
B(2, 1:4) = [-cline_para.alpha3D-2*cline_para.beta3D, 2*cline_para.alpha3D+5*cline_para.beta3D, -cline_para.alpha3D-4*cline_para.beta3D, cline_para.beta3D];
for i=3:stack_z_size-2
    B(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end
B(stack_z_size-1, stack_z_size-3:stack_z_size) = ...
    [cline_para.beta3D, -cline_para.alpha3D-4*cline_para.beta3D, 2*cline_para.alpha3D+5*cline_para.beta3D, -cline_para.alpha3D-2*cline_para.beta3D];
B(stack_z_size, stack_z_size-2:stack_z_size) = ...
    [cline_para.beta3D, -2*cline_para.beta3D, cline_para.beta3D];
[L, U] = lu(B + cline_para.gamma .* eye(stack_z_size,stack_z_size));
Binv = inv(U) * inv(L); % Computing Ainv using LU factorization




        %% gradient
        
        [fy, fx] = gradient(im*cline_para.gradient_force); %computing the gradient
      %  fx=zeros(size(fx));
%         fy=fy;
if cline_para.show_flag
         figure;
         for ishow=1:4
         subplot(2,2,ishow)
 imagesc(im(:,:,round(1+(ishow-1)*(stack_z_size-1)/3)));
    hold on
   
         end
 
end
        %% prepare for cline_para.iterations
        xyzs_avg = zeros(size(xyzs));
        counter = 0;
        %moving the snake in each iteration
        [~,~,tgrid]=meshgrid(1,1:size(xyzs,1),1:stack_z_size);
        for i=1:cline_para.iterations;
            %calculate image gradient force on the nodes
             % actual interpolation, slow!!!
                xs=(xyzs(:,1,:));
                ys=(xyzs(:,2,:));
                
                fs = cat(2,interp3(fy,xs,ys,tgrid,'*linear'),...
                    interp3(fx,xs,ys,tgrid,'*linear'))    ;     
                
%                 
                 fs(isnan(fs)) = 0;
%                 if i<cline_para.iterations/2
%                 fs(end,1)=(xs(end-5)-xs(end))*.3;
%                 fs(1,1)=(xs(5)-xs(1))*.3;
%                 end
% 
% %  fs(end,:)=1*nanmean(fs(round(length(fs)/3):end,:));
%    fs(1,:)=1*nanmean(fs(1:round(length(fs)/3),:));

            %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
            % viscos force between frames on all nodes:        
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
            
timeForce=zeros(size(xyzs ));
%timeForce(:,:,2:end-1)=diff(xyzs,2,3);

%X relax
            
             tempOut= Ainv * squeeze(cline_para.gamma*xyzs(:,1,:)+ cline_para.kappa*fs(:,1,:)...
                 +cline_para.alpha3D*timeForce(:,1,:));
            
           xyzs(:,1,:) = reshape(tempOut,size(xyzs(:,1,:)));
            %yRelax
%                tempOut= Ainv * squeeze(cline_para.gamma*xyzs(:,2,:) + cline_para.kappa*fs(:,2,:)...
%                 +cline_para.alpha3D*timeForce(:,2,:));
%                         xyzs(:,2,:) = reshape(tempOut,size(xyzs(:,1,:)));
        %timeRelax
           tempOut= (Binv * squeeze(cline_para.gamma*xyzs(:,1,:))')';
            
            xyzs(:,1,:) = reshape(tempOut,size(xyzs(:,1,:)));
        
            
            % take average value at the end of the cline_para.iterations
            if i > cline_para.iterations - 20
                xyzs_avg = xyzs_avg + xyzs;
                counter = counter + 1;
            end
            
            %end points are fixed

   %          xyzs(1,2,:)=1; 
  %           xyzs(end,2,:)=row;
            L=size(xyzs,1)-1;
            %redistribute points evenly
 s=cat(1,zeros(1,1,stack_z_size),cumsum(sqrt(diff(xyzs(:,1,:),[],1).^2+diff(xyzs(:,2,:),[],1).^2),1));
  s=squeeze(s);
   lgrid=(0:1/L:1)'*squeeze(max(s,[],1));
   

lcell=mat2cell(squeeze(lgrid),row,ones(1,stack_z_size));
scell=mat2cell(s,row,ones(1,stack_z_size));
xcell=mat2cell(squeeze(xyzs(:,1,:)),row,ones(1,stack_z_size));
ycell=mat2cell(squeeze(xyzs(:,2,:)),row,ones(1,stack_z_size));

xcell=cell2mat(cellfun(@interp1,scell,xcell,lcell,'uniformoutput',0));
ycell=cell2mat(cellfun(@interp1,scell,ycell,lcell,'uniformoutput',0));

 xyzs(:,1,:)=reshape(xcell,row,1,stack_z_size);
    xyzs(:,2,:)=reshape(ycell,row,1,stack_z_size);
%           
    %view contour
    if cline_para.show_flag
        
    if exist('h','var')
        for ishow=1:4
    set(h(ishow),'Xdata',xyzs(:,1,round(1+(ishow-1)*(stack_z_size-1)/3)),...
        'Ydata' ,xyzs(:,2,round(1+(ishow-1)*(stack_z_size-1)/3)))
        end
    else
        
         for ishow=1:4
         subplot(2,2,ishow)
         h(ishow)=scatter(xyzs(:,1,round(1+(ishow-1)*(stack_z_size-1)/3)),...
             xyzs(:,2,round(1+(ishow-1)*(stack_z_size-1)/3)));
         
         end
    
    end
    end
    drawnow
%     
          
        end
         %avg last% 50 iterations
         xyzs_avg = xyzs_avg / counter;

end
