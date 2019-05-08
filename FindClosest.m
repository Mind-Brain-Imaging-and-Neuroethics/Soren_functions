function [indicesout,itemsout] = FindClosest(list,items,first)

if nargin < 3
    first = 1; 
end

for c = 1:length(items)
   tmp = find(abs(list-items(c)) == min(abs(list-items(c))));
    indicesout(c) = tmp(first);
end

itemsout = list(indicesout);