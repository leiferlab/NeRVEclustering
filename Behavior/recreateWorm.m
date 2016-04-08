function centerline_out=recreateWorm(varargin)
% recreate Worm takes wormcentered data several coordinates of point
% indices in order to recreate the shape of the worm. 

%%% centerline_out=recreateWorm(wc,ref_index,ref_point,ref_angle,worm_length)
%%% wc : worm centered coordinates
%%% ref_point : coordinate of reference point used to position worm
%%% ref_angle :angle of vector at that point
%%% ref_index : index of point being used as the refence
%%% worm_length : total length of the worm

%%% centerline_out=recreateWorm(wc,ref_index,ref_point)
%%% wc : worm centered coordinates
%%% ref_point : coordinate of reference point used to position worm (2 of
%%% them)
%%% ref_index : index of point being used as the refence (2 of them)

wc=varargin{1};
ref_point=varargin{3};
ref_index=varargin{2};

if nargin==5
    ref_angle=varargin{4};
worm_length=varargin{end};

angleOffset=ref_angle-wc(ref_index);
wc2=wc+angleOffset;

tVector=worm_length/length(wc2)*[sin(wc2) cos(wc2)];
X=cumsum(tVector);
centerline_out=bsxfun(@plus, X, -X(ref_index,:)+ref_point);
elseif nargin>=3

        
        
    delta_ref=ref_point(2,:)-ref_point(1,:);
    ref_angle=atan2(delta_ref(2),delta_ref(1));
    tVector=[sin(wc(:)) cos(wc(:))];
X=cumsum(tVector);
delta_x_ref=X(ref_index(2),:)-X(ref_index(1),:);
sample_angle=atan2(delta_x_ref(2),delta_x_ref(1));


angleOffset=ref_angle-sample_angle;

wc2=wc-angleOffset;

tVector=[sin(wc2(:)) cos(wc2(:))];
X=cumsum(tVector);
        if nargin==4;
            worm_length=varargin{end};
            X=X*worm_length/length(wc2);
        else
            delta_x_ref=X(ref_index(2),:)-X(ref_index(1),:);
X=X*mean(delta_ref./delta_x_ref);

        end
        

centerline_out=bsxfun(@plus, X, -X(ref_index(1),:)+ref_point(1,:));
end

    
    