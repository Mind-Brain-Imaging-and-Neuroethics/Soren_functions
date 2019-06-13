function [newlim] = Set_Clim(fighandle,clim)

ax = findall(fighandle,'Type','Axes');

for c = 1:length(ax)
    try
        Clims(c,:) = get(ax(c),'CLim');
    end
end

for c = 1:length(ax)
    try
        set(ax,'CLim',clim)
    end
end