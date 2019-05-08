function [fieldout] = getfield_nest(structin,field)

newfield = strsplit(field,'.');

streval = [];
for c = 1:length(newfield)
streval = [streval '.(''' newfield{c} ''')'];
end

eval(['fieldout = structin' streval ';']);