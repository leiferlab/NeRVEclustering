function [FX,FY,FZ]=fiducials2mat(fiducialPoints,n)

select=cellfun(@(x) ~isempty(x),fiducialPoints);
temp=cellfun(@(x) (x(:,1:3)),fiducialPoints(select),'uniform',0);
temp=[temp{:}];
empty=cellfun(@(x) isempty(x),temp);
temp(empty)={nan};
temp=cell2mat(temp);

FX=nan(size(temp,1),length(fiducialPoints));
FY=nan(size(temp,1),length(fiducialPoints));
FZ=nan(size(temp,1),length(fiducialPoints));
FX(:,select)=temp(:,1:3:end);
FY(:,select)=temp(:,2:3:end);
FZ(:,select)=temp(:,3:3:end);
if nargin>1
    FX=FX(1:n,:);
    FY=FY(1:n,:);
    FZ=FZ(1:n,:);
end
