function [corrrho,corrp,h,f,tb] = nicecorrplot(a,b,labels,varargin)

if EasyParse(varargin,'RemoveOutliers','on')
    rmoutliers = [find(isoutlier(a))' find(isoutlier(b))'];
    
    a(rmoutliers) = [];
    b(rmoutliers) = [];
end

if size(a,1) == 1
    a = a';
end

if size(b,1) == 1
    b = b';
end

[B,~,~,~,stats] = regress(b,horzcat(ones(length(a),1),a));

if EasyParse(varargin,'Cov')
[corrrho,corrp] = partialcorr(a,b,EasyParse(varargin,'Cov'),'Type','Spearman')
else
[corrrho,corrp] = corr(a,b,'Type','Spearman')
end

h = scatter(a,b,64,'filled');
xl = xlabel(labels{1},'FontSize',16);
yl = ylabel(labels{2},'FontSize',16);
hold on;
f = plot(linspace(min(a),max(a),1000),B(1)+B(2)*linspace(min(a),max(a),1000),'r--');
pos = get(gca,'Position');
if CheckInput(varargin,'Plot') && EasyParse(varargin,'Plot','r')
tb = annotation('textbox',[pos(1)+0.2 pos(2)+0.4 0.3 0.3],'String',{['rho = ' num2str(round(corrrho,3))]},'FitBoxToText','on');
else
tb = annotation('textbox',[pos(1)+0.2 pos(2)+0.5 0.3 0.3],'String',{['rho = ' num2str(round(corrrho,3))],['p = ' num2str(round(corrp,3,'significant'))]},'FitBoxToText','on');
end
%ylim([-0.1 0.1])

set(f,'LineWidth',1.5)
set(tb,'LineStyle','none','FontSize',16)

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'fontsize',14,...
  'LineWidth'   , 1         );

set(xl,'FontSize',20)
set(yl,'FontSize',20)
