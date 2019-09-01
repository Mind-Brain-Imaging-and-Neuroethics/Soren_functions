function [newlim] = Set_Clim(fighandle,ylim)

ax = findall(fighandle,'Type','Axes');

% for c = 1:length(ax)
%     try
%         Ylims(c,:) = get(ax(c),'CLim');
%     end
% end

for c = 1:length(ax)
    try
        set(ax,'YLim',ylim)
    end
end