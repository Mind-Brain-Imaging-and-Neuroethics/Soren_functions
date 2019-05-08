function stdshade(F,amatrix,acolor,alpha,dimn,method,smth)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data.
%sem/std is shown as shading.
% - dimn specifies the dimension over which to plot the data (default = 2 -
% rows are observations)
% - acolor defines the used color (default is red)
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

%modified by SWT on 19/03/2019

if dimn == 1
    amatrix = amatrix';
end

if ~exist('method','var')
    method = 'sem';
end

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r';
end

if exist('F','var')==0 || isempty(F);
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;
end

if ne(size(F,1),1)
    F=F';
end

amean=smooth(nanmean(amatrix),smth)';

switch method
    case 'std'
        astd=nanstd(amatrix); % to get std shading
    case 'sem'
        astd=nanstd(amatrix)/sqrt(size(amatrix,1)-1); % to get sem shading
    case 'mad'
        astd = mad(amatrix);
end

if exist('alpha','var')==0 || isempty(alpha)
    fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none','HandleVisibility','off');
    acolor='k';
else fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none','HandleVisibility','off');
end

if ishold==0
    check=true; else check=false;
end

hold on;
if ischar(acolor)
    plot(F,amean,acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
else
    plot(F,amean,'Color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
end


if check
    hold off;
end

end



