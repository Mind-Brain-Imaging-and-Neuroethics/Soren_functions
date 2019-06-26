function merged = mergestructs(structsin)

if ~iscell(structsin)
    fields = cell_unpack(fieldnames_recurse(structsin(1)));
    merged = struct;
    for c = 1:length(fields)
        dimn = size(getfield_nest(structsin(1),fields{c}));
        if ismember(1,dimn)
            dimn = find(dimn == 1,1);
        else
            dimn = length(dimn)+1;
        end
        alldata = cell(1,length(structsin));
        for cc = 1:length(structsin)
            alldata{cc} = getfield_nest(structsin(cc),fields{c});
        end
        merged = assignfield_nest(merged,fields{c},cat(dimn,alldata{:}));
    end
else
    fields = cell_unpack(fieldnames_recurse(structsin{1}));
    merged = struct;
    for c = 1:length(fields)
        dimn = size(getfield_nest(structsin{1},fields{c}));
        if ismember(1,dimn)
            dimn = find(dimn == 1,1);
        else
            dimn = length(dimn)+1;
        end
        alldata = cell(1,length(structsin));
        for cc = 1:length(structsin)
            alldata{cc} = getfield_nest(structsin{cc},fields{c});
        end
        merged = assignfield_nest(merged,fields{c},cat(dimn,alldata{:}));
    end
end
