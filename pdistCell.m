function DmatAll=pdistCell(pointsCell,pointsIdx,dims)
% takes cell array of coordinates and makes a square pdistance for each of
% them. Only uses points specified by pointsIdx and dims

%default use 2d
if nargin<=2
    dims=2;
end
if nargin==1
    pointsIdx=[];
end

%default pointsIdx is all points in first entry of pointsCell
if length(dims)==1
    dims=1:dims;
end
if isempty(pointsIdx)
  pointsIdx=1:size(  pointsCell{1},1);
end
DmatAll=nan(length(pointsIdx),length(pointsIdx),length(pointsCell));

for i=1:length(pointsCell);
    points=pointsCell{i};
    if ~isempty(points)
    points=points(:,dims);
    
    emptyX=cellfun(@(x) ~isempty(x),points(:,1));
plotIdx=find(emptyX);
%if ii==1
X0=nan(length(pointsIdx));
%end

  [ plotIdx2,ia,ic]= intersect(plotIdx,pointsIdx);
   points=points(plotIdx2,:);
   points=cell2mat(points);
   
    Dmat=squareform(pdist(points));
    end
    if ~isempty(Dmat);
X0(ic,ic)=Dmat; 
    end
    DmatAll(:,:,i)=X0;
   
    
    
end
