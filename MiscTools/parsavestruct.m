function parsavestruct(fname, outputData)
fieldstring=fieldnames(outputData);
for i=1:length(fieldstring)
    eval([fieldstring{i} '= outputData.' fieldstring{i} ';']);

    
    try
        save(fname,fieldstring{i},'-append');
        
    catch
        save(fname,fieldstring{i});
    end
    
end



