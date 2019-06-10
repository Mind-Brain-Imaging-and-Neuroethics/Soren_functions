function AddFigureLabel(ax,string,topo)

outpos = ax.OuterPosition;
inpos = ax.Position;

a = annotation(gcf,'textbox',[0 0 0.02 0.075],...
    'String',string,...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'FontSize',24,'FontWeight','bold'); %just get the textbox size from this, then clear it

tbsize = get(a,'Position');
tbsize = tbsize(3:4);
delete(a)


if ~exist('topo','var') || strcmpi(topo,'no')
offset = [outpos(1) inpos(2)+inpos(4)-tbsize(2)+1/2*(outpos(2)+outpos(4)-inpos(2)-inpos(4))];
else
    offset = [outpos(1)-tbsize(1)*2 inpos(2)+inpos(4)+-tbsize(2)/2+1/2*(outpos(2)+outpos(4)-inpos(2)-inpos(4))];
end

offset(find(offset < 0)) = 0;

a = annotation(gcf,'textbox',[offset 0.02 0.075],...
    'String',string,...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'FontSize',24,'FontWeight','bold'); %just get the textbox size from this, then clear it





