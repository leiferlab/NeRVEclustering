killIdx=[1 5];

for i=1:length(fiducialPoints)
        temp=fiducialPoints{i};

    for iKill=killIdx
    temp(iKill,:)={[],[],[],[],[]};
    end
        fiducialPoints{i}=temp;

end
