function Plot_sigmask(ax,sigmask,type,varargin)

hold(gca,'on')

if nargin < 3
    type = 'bar';
end

if ~isempty(find(size(sigmask) == 1))
   sigmask = horz(sigmask); 
end

%sigmask = sigmask > 0; % just to make sure it's zeros and ones

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
varargin = setdefault(varargin,'cmap',gray);
varargin = setdefault(varargin,'linewidth',2);

color = EasyParse(varargin,'color');
cmap = EasyParse(varargin,'cmap');
alpha = EasyParse(varargin,'alpha');
linewidth = EasyParse(varargin,'linewidth');

xlim = get(gca,'XLim');
xax = linspace(xlim(1),xlim(2),size(sigmask,2)); %time points must be the second dimension of sigmask

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
    if length(stopindices) < length(startindices)
       stopindices = [stopindices length(sigmask)]; 
    end
    indices = [vert(startindices) vert(stopindices)];
    
    if ~isempty(indices)
    for c = 1:size(indices,1)
        line([xax(indices(c,1)) xax(indices(c,2))],[yl(2)-ysize*0.04 yl(2)-ysize*0.04],...
            'color',[color alpha],'linewidth',linewidth);
    end
    end
elseif strcmpi(type,'cmapline') % for plotting the number of significant channels as the colour of the line
    yl = get(ax,'YLim');
    ysize = yl(2)-yl(1);
    y = repmat([yl(2)-ysize*0.04],size(xax));
    z = zeros(size(xax));
    col = sum(sigmask,find(size(sigmask) ~= length(xax)));
    cindex = linspace(size(sigmask,1),0,64);
    col = cmap(FindClosest(cindex,col,1),:);
    col = permute(col,[3 1 2]);
    
    surface([xax;xax],[y;y],[z;z],[col;col],...
        'facecolor','no','edgecolor','interp','linewidth',linewidth)
    %colormap(cmap)
end