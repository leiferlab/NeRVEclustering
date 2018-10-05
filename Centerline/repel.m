function repforce = repel(P,cd, mu)

nPoints = length(P);

L = P(1:ceil(nPoints/2), :);
M = P(nPoints - floor(nPoints/2):nPoints, :);


repel1 = zeros(length(L), 2);
force = zeros(nPoints,2);
for i = 1:length(L(:,1))
    FL = [0 0];
    for j = 1:length(M)
        fj = sqrt(sum(((L(i,:)-M(j,:)).^2)))/cd;
        pickone = max([0, 1-fj]);
        fx = pickone*((L(i,1)-M(j,1))/cd);
        fy = pickone*((L(i,2)-M(j,2))/cd);
        fxy = [fx fy];
        
        
        FL = FL + fxy;
    end
    
    repel1(i,:) = FL;
  
    
    
end
repel2 = zeros(length(L), 2);
for i = 1:length(M)
    FL = [0 0];
    for j = 1:length(L(:,1))
        fj = sqrt(sum(((M(i, :)-L(j,:)).^2)))/cd;
        pickone = max([0, 1-fj]);
        fx = pickone*((M(i,1)-L(j,1))/cd);
        fy = pickone*((M(i,2)-L(j,2))/cd);
        fxy = [fx fy];
        
        
        FL = FL + fxy;
    end
    
    repel2(i,:) = FL;
  
    
    
end

force(1:ceil(nPoints/2),:) = repel1;
force(nPoints-50:nPoints,:) = repel2(end-50:end,:);

repforce = mu*force;
end