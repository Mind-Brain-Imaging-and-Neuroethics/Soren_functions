function Plot_sigmask(ax,sigmask,type,varargin)

hold(gca,'on')

if nargin < 3
    type = 'bar';
end

sigmask = sigmask > 0; % just to make sure it's zeros and ones

if strcmpi(type,'bar')
    varargin = setdefault(varargin,'color',[0.5 0.5 0.5]);
else
    varargin = setdefault(varargin,'color',[0.15 0.15 0.15]);
end
if strcmpi(type,'bar')
    varargin = setdefault(varargin,'alpha',0.5);
else
    varargin = setdefault(varargin,'alpha',0.10);
end
varargin = setdefault(varargin,'linewidth',2);

color = EasyParse(varargin,'color');
alpha = EasyParse(varargin,'alpha');
linewidth = EasyParse(varargin,'linewidth');

xlim = get(gca,'XLim');
xax = linspace(xlim(1),xlim(2),length(sigmask));

if strcmpi(type,'shade')
    yl = get(ax,'YLim');
    patchstep = sigmask*(yl(2)-yl(1));
    patchstep = patchstep + yl(1);
    area(xax,patchstep,yl(1),'FaceAlpha',alpha,'LineStyle','none','FaceColor',color,'HandleVisibility','off')
    set(ax,'YLim',yl)
elseif strcmpi(type,'bar')
    yl = get(ax,'YLim');
    ysize = yl(2)-yl(1);
    diffmask = diff([0 sigmask]);
    startindices = find(diffmask == 1);
    stopindices = find(diffmask == -1);
    indices = [vert(startindices) vert(stopindices)];
    
    for c = 1:size(indices,1)
        line([xax(indices(c,1)) xax(indices(c,2))],[yl(2)-ysize*0.04 yl(2)-ysize*0.04],...
            'color',[color alpha],'linewidth',linewidth);
    end
end