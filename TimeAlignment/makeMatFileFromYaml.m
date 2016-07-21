function makeMatFileFromYaml(yamlFile)

if nargin==0
    yamlFile=uipickfiles('filterspec','E:\');
end
for iYaml=1:length(yamlFile)
mcdf=Mcd_Frame;
mcdf=mcdf.yaml2matlab(yamlFile{iYaml});
yamlOut=strrep(yamlFile{iYaml},'.yaml','YAML');
save(yamlOut,'mcdf');
clear mcdf
end
