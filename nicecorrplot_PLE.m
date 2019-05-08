function [h,f,tb] = nicecorrplot(a,b)

figure

rmoutliers = [find(isoutlier(a))' find(isoutlier(b))'];

a(rmoutliers) = [];
b(rmoutliers) = [];


[B,~,~,~,stats] = regress(b,horzcat(ones(length(a),1),a));

[corrrho,corrp] = corr(a,b,'Type','Spearman')

h = scatter(a,b,'filled')
xlabel('PLE','FontSize',16)
ylabel('Slope of self effect over delays','FontSize',16)
hold on;
f = plot(linspace(min(a),max(a),1000),B(1)+B(2)*linspace(min(a),max(a),1000),'r--')
tb = annotation('textbox',[0.2 0.6 0.3 0.3],'String',{['rho = ' num2str(round(corrrho,3))],['p = ' num2str(round(corrp,3,'significant')) '*']},'FitBoxToText','on');
ylim([-0.1 0.1])

set(f,'LineWidth',1.5)
set(tb,'LineStyle','none','FontSize',14)

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'XTick',0.7:0.2:1.5,...
  'YTick',-0.1:0.04:0.1,...
  'LineWidth'   , 1         );