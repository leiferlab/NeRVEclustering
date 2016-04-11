function [centerline_out,err,wc]=EigenwormApprox(centerline_in,refIdx2,eigbasis,eigen_range)
%% recreates a centerline using the first 6 eigen worms
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
