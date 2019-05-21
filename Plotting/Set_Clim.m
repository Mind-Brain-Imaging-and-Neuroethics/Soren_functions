function [newlim] = Set_Clim(fighandle,clim)

ax = findall(fighandle,'Type','Axes');

for c = 1:length(ax)
   Clims(c,:) = get(ax(c),'CLim');
end

for c = 1:length(ax)
   set(ax,'CLim',clim) 
end