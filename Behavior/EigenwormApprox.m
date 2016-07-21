function [centerline_out,err,wc]=EigenwormApprox(centerline_in,refIdx2,eigbasis,eigen_range)
% recreates a centerline using the first 6 eigen worms. The 2 points of the
% intial centerline, given by refIdx2, are pinned to constrain the fitting.
% 

%centerline_in : the initial centerline to be reapproximated
%refIdx2: the points of the initial worm to fix
% eigenbasis: the matrix of eigenworms to be used to reconstruct the
% centerline in, as nWorms x 100
% eigen_range: the range of centerline points to use as the entire length
% of the worm. ex. 1:50 will say that point 50 is the end of the worm.

if nargin==3
    eigen_range=1:(length(centerline_in)-1);
end
%centerline_in=distanceInterp(centerline_in,length(centerline_in));

wc=FindWormCentered(centerline_in);
weights=sqrt(sum(eigbasis(:,eigen_range).^2,2));
wc2=((eigbasis(:,eigen_range)*wc(eigen_range,:))./weights)'*eigbasis(:,1:end-1);
ref_points=centerline_in(refIdx2,:);
centerline_out=recreateWorm(wc2,refIdx2,ref_points);
%centerline_out=distanceInterp(centerline_out,length(centerline_in));
if nargout>1
err=sum(sqrt(sum((centerline_in(eigen_range,:)-centerline_out(eigen_range,:)).^2,2)))/length(eigen_range);
end
% plot(centerline_in(:,1),centerline_in(:,2));
% hold on
% plot(centerline_out(:,1),centerline_out(:,2));
% hold off
