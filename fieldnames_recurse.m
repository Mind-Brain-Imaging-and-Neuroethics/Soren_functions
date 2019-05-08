function [fields] = fieldnames_recurse(structin)

fields = fieldnames(structin);

for c = 1:length(fields)
    if isstruct(structin.(fields{c}))
        tmp = fieldnames_recurse(structin.(fields{c}));
        fields{c} = cellcat(fields{c},tmp,'.');
    end
end