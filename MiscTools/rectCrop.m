function imOut=rectCrop(imIn,rect1,padVal)
if nargin==2
    padVal=0;
end


if all(rect1>=0) && all(size(imIn)>=rect1([4,3]));
imOut=imIn((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
else
    if any(rect1<0);
        padSize=-min(rect1);
    else padSize=0;
    end
    
    if any(rect1([4,3])>=size(imIn))
    padSize=max([rect1([4,3])-size(imIn),padSize]);
    end
    
    imIn=padarray(imIn,[padSize,padSize],padVal,'both');
    rect1=rect1+padSize;
imOut=imIn((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
    
    
end

