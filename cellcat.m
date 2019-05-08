function [cellout] = cellcat(stringin,cellin,joiner)
    for c = 1:length(cellin)
        if iscell(cellin{c})
            cellout{c} = cellcat(stringin,cellin{c},joiner);
        else
        cellout{c} = [stringin joiner cellin{c}];
        end
    end
end