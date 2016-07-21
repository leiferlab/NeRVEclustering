function parsave(fname, outputData,outputString)
eval([outputString '= outputData;']);
try
    
save(fname,outputString,'-append');
catch
    save(fname,outputString);
end

end

