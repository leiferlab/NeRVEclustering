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

%get inputs
wc=varargin{1};
ref_point=varargin{3};
ref_index=varargin{2};

if nargin==5
    ref_angle=varargin{4};
    worm_length=varargin{end};
    %get angle offset by finding difference between the ref_angle and the
    %wormcentered angle at that point
    angleOffset=ref_angle-wc(ref_index);
    %add that offset everywhere
    wc2=wc+angleOffset;
    %build all tangent vectors using trig, also rescale length
    tVector=worm_length/length(wc2)*[sin(wc2) cos(wc2)];
    %add up all the steps to build up centerline coordinates
    centerline=cumsum(tVector);
elseif nargin>=3
    %find difference between two reference points
    delta_ref=ref_point(2,:)-ref_point(1,:);
    %get angle betwen the two ref points
    ref_angle=atan2(delta_ref(2),delta_ref(1));
    %get all vectors for the centerline based on worm centered coordinate
    tVector=[sin(wc(:)) cos(wc(:))];
    
    centerline_raw=cumsum(tVector);
    %find angle between ref ponits without scaling/rotation yet
    delta_x_ref=centerline_raw(ref_index(2),:)-centerline_raw(ref_index(1),:);
    sample_angle=atan2(delta_x_ref(2),delta_x_ref(1));
    %find difference between angles in centerline and the ref points, will need
    %to rotate the centerline to match that
    angleOffset=ref_angle-sample_angle;
    %rotate centerline by adding agle offset to wormcentered
    wc2=wc-angleOffset;
    %build centerline again now with proper angle
    tVector=[sin(wc2(:)) cos(wc2(:))];
    centerline=cumsum(tVector);
    if nargin==4; %scale centerline by length if total length given
        worm_length=varargin{end};
        centerline=centerline*worm_length/length(wc2);
    else
        %if total length not given, scale up by difference between ref poitns
        delta_x_ref=centerline(ref_index(2),:)-centerline(ref_index(1),:);
        centerline=centerline*mean(delta_ref./delta_x_ref);
    end
end

%add xy offset 
xy_offset= -centerline(ref_index(1),:)+ref_point(1,:);
centerline_out=bsxfun(@plus, centerline, xy_offset);

