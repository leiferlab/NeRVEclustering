function P=WormBrain3Points(WormProj)

WormProj=pedistalSubtract(WormProj);
WormProj=smooth2a(WormProj,7,7);
WormProj=normalizeRange(WormProj);

for iterations=1:2
botWormProjBW=im2bw(WormProj,(1-.1*iterations)*graythresh(WormProj));
botWormProjBW=imopen(botWormProjBW,true(5));

botLine=botWormProjBW(end,:);
botWormProjBW=watershedFilter(botWormProjBW,5);
botWormProjBW=imclearborder(botWormProjBW,4);

botStats=regionprops(botWormProjBW,'Area','Centroid');
blobAreas=[botStats.Area]';
centroids=cell2mat({botStats.Centroid}');
%remove centroids more than 70% across the image

if ~isempty(centroids);
select=centroids(:,1)<size(WormProj,2)*.7;
centroids=centroids(select,:);
blobAreas=blobAreas(select);
end

if iterations>1 && isempty(centroids);
    P=P(all(P,2),:);
    return
end
cPoint=mean(find(botLine));
if isnan(cPoint);
    cPoint=size(WormProj,2)/2;
end

P(3,:)=[cPoint,size(WormProj,1)];

if ~isempty(centroids)
P(2,:)=centroids(blobAreas==max(blobAreas),:);
posDis=sqrt(sum(bsxfun(@minus,centroids,P(3,:)).^2,2));
P(1,:)=centroids(posDis==max(posDis),:);
else
    P(2,:)=[cPoint,size(WormProj,1)*2/3];
end


if all(P(1,:)==P(2,:))
P(1,:)=[];
end

if length(P)>2
        P=P(all(P,2),:);

    return
end

end


