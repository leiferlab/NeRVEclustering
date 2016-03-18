

function fiducialPoints=deleteNeuron(fiducialPoints,killIdx);

for i=1:length(fiducialPoints)
        temp=fiducialPoints{i};
if ~isempty(temp)
    for iKill=killIdx
    temp(iKill,:)=cell(1,size(temp,2));
    end
end
        fiducialPoints{i}=temp;

end
