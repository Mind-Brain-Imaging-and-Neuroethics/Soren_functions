function [newlim] = Normalize_Clim(fighandle,equal)

ax = findall(fighandle,'Type','Axes');

for c = 1:length(ax)
   Clims(c,:) = get(ax(c),'CLim');
end

newlim = [min(Clims(:,1)) max(Clims(:,2))];

if nargin > 2 || equal
    maxabs = max(abs(newlim));
    newlim = [-maxabs maxabs];
end
for c = 1:length(ax)
   set(ax,'CLim',newlim) 
end