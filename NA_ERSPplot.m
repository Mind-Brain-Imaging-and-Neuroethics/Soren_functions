function NA_ERSPplot(ersp,orig_trange,trange,pdifs,varargin)

xax = linspace(1,orig_trange(2)-orig_trange(1),length(trange));

plot(xax,mean(mean(ersp{1},3),1),'b')

hold on

plot(xax,mean(mean(ersp{2},3),1),'r')

if CheckInput(varargin,'PlotPseudo')
    pseudoersp = EasyParse(varargin,'PlotPseudo')
   plot(xax,mean(mean(pseudoersp{1}(:,1:length(trange),:),3),1),'b--')
   plot(xax,mean(mean(pseudoersp{2}(:,1:length(trange),:),3),1),'r--')
end

patchindex = pdifs < 0.05;
yl = ylim;
patchstep = patchindex*(yl(2)-yl(1));
patchstep = patchstep + yl(1);
area(xax,patchstep,yl(1),'FaceAlpha',0.3,'LineStyle','none','FaceColor',[0.5 0.5 0.5])
ylim(yl)
%legend({'Corrected Prestim Low','Corrected Prestim High','Significant Differences'})
xlabel('Time (ms)')
ylabel('ERSP')