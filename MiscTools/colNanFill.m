function stackIn=colNanFill(stackIn)
% colNanFill takes an image stack with nans due to warping and fills them
% in with the nears non nan value in that column.

 imSize=size(stackIn);
 nanMap=isnan(stackIn);

for i=1:imSize(2)
    for j=1:imSize(3)
        tempCol=stackIn(:,i,j);
        
        nanTemp=nanMap(:,i,j);
        tempCol(nanTemp)=0;
        firstVal=find(tempCol,1,'first');
        lastVal=find(tempCol,1,'last');
        
        tempCol(1:firstVal)=tempCol(firstVal);
        tempCol(lastVal:end)=tempCol(lastVal);
        stackIn(:,i,j)=tempCol;
    end
end







