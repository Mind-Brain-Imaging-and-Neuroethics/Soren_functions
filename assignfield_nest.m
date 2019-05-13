function [structin] = assignfield_nest(structin,field,data)

newfield = strsplit(field,'.');

streval = [];
for c = 1:length(newfield)
streval = [streval '.(''' newfield{c} ''')'];
end

eval(['structin' streval ' = data;']);