function AddFigureLabel(ax,string,postype,varargin)

outpos = ax.OuterPosition;
inpos = ax.Position;

if nargin < 3
   postype = 'upperleft'; 
end

varargin = setdefault(varargin,'FontSize',24);

a = annotation(gcf,'textbox',[0 0 0.02 0.075],...
    'String',string,...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'FontSize',EasyParse(varargin,'FontSize'),'FontWeight','bold'); %just get the textbox size from this, then clear it

tbsize = get(a,'Position');
tbsize = tbsize(3:4);
delete(a)


if strcmpi(postype,'topo_upperleft')
    offset = [outpos(1) inpos(2)+inpos(4)-tbsize(2)+1/2*(outpos(2)+outpos(4)-inpos(2)-inpos(4))];
elseif strcmpi(postype,'upperleft')
    offset = [outpos(1)-tbsize(1)*2 inpos(2)+inpos(4)+-tbsize(2)/2+1/2*(outpos(2)+outpos(4)-inpos(2)-inpos(4))];
elseif strcmpi(postype,'middle_left')
    offset = [outpos(1)-tbsize(1) inpos(2)+0.5*inpos(4)+-tbsize(2)/2];
end

offset(find(offset < 0)) = 0;

a = annotation(gcf,'textbox',[offset 0.02 0.075],...
    'String',string,...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'FontSize',EasyParse(varargin,'FontSize'),'FontWeight','bold'); 





