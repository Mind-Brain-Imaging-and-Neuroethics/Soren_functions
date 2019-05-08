function [boundaries] = GetGroupBoundaries(groupSeries)
%returns the boundaries from a sequence
%each sequence must be a COLUMN of the groupSeries matrix

len = size(groupSeries,1);

boundaries = zeros(len-1,size(groupSeries,2));

for q = 1:size(groupSeries,2)
    for c = 1:len-1
        if groupSeries(c+1,q) ~= groupSeries(c,q) && groupSeries(c,q) ~= 0
            boundaries(c,q) = 1;
        end
    end
end
%boundaries = squeeze(boundaries);