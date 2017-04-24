function wormcentered=FindWormCentered(centerline)

%FindWormCentered takes centereline coordinates (100x2x time) and
%calculates the wormcentered coordinatesystem, from A Linder. 
%%
numframes = 1:length(centerline(1,1,:));


numcurvpts=size(centerline,1)-2; 

tanvecs = zeros(numcurvpts+1, 2, length(numframes));
atdf2 = zeros(numcurvpts+1, length(numframes));
deltaTheta = zeros(numcurvpts, length(numframes));

wormcentered = zeros(numcurvpts+1, length(numframes));
for j = 1:length(numframes)
    df2 = diff(centerline(:,:,j),1,1); % calculate tangent vector along curve (not normalized)
    
    tanvecs(:,:, j) = df2;
    
    angles = unwrap(atan2(-df2(:,2), df2(:,1))); % angle of tangent vector.  unwrapped to avoid 2pi jumps
    
    
     atdf2(:,j) = angles;
    %deltaTheta(:,j) = unwrap(diff(atdf2,1)); % derivative of angle with respect to path length
 
    %subtract average theta to get 'wormcentered' coords
    avg = mean(angles);
    
    wormcentered(:,j) = angles - avg;
end