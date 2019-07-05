function compstats = NA_compare_outputs(settings,allmeas1,allmeas2)

opts.nrand = 1000;
opts.minnbchan = 1;
opts.parpool = settings.pool;
opts.eqinterval = [-2 2];

for c = 1:settings.nfreqs
    %pt_diff_stats{c} = EasyClusterCorrect({permute(squeeze(allmeas1{c}.naddersp.diff(:,:,2,:)-allmeas1{c}.naddersp.diff(:,:,1,:)),[1 3 2]),...
    %    permute(squeeze(allmeas2{c}.naddersp.diff(:,:,2,:)-allmeas2{c}.naddersp.diff(:,:,1,:)),[1 3 2])},...
    %    settings.datasetinfo,'ft_statfun_fast_signrank',opts);
    %ttv_diff_stats{c} = EasyClusterCorrect({permute(allmeas1{c}.ttversp.real,[1 3 2]),permute(allmeas2{c}.ttversp.real,[1 3 2])},...
    %    settings.datasetinfo,'ft_statfun_fast_signrank',opts);

    pt_diff_stats{c} = EasyClusterCorrect({permute(squeeze(allmeas1{c}.naddersp.diff(:,:,2,:)-allmeas1{c}.naddersp.diff(:,:,1,:)),[1 3 2]),...
        permute(squeeze(allmeas2{c}.naddersp.diff(:,:,2,:)-allmeas2{c}.naddersp.diff(:,:,1,:)),[1 3 2])},...
        settings.datasetinfo,'ft_statfun_tost',opts);
    ttv_diff_stats{c} = EasyClusterCorrect({permute(allmeas1{c}.ttversp.real,[1 3 2]),permute(allmeas2{c}.ttversp.real,[1 3 2])},...
        settings.datasetinfo,'ft_statfun_tost',opts);        


    ersp_diff_stats{c} = EasyClusterCorrect({permute(allmeas1{c}.ersp.real,[1 3 2]) permute(allmeas2{c}.ersp.real,[1 3 2])},...
        settings.datasetinfo,'ft_statfun_fast_signrank',opts);
end

compstats.ptdiff = pt_diff_stats;
compstats.ttvdiff = ttv_diff_stats;
compstats.erspdiff = ersp_diff_stats;

save(fullfile(settings.outputdir,[settings.datasetname '_compstats.mat']),'compstats','-v7.3')

cd(settings.outputdir)

figure

p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*4 pos(4)*4],'Color','w');

p.pack('v',{1/2 1/2})

p(1).pack('h',repmat({1/settings.nfreqs},1,settings.nfreqs))

p(2).pack('h',repmat({1/settings.nfreqs},1,settings.nfreqs))


for c = 1:settings.nfreqs
    
    %Pseudotrial-based differences
    
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    p(1,c).pack();
    for cc = 1:4
        p(1,c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(1,c,1).select();
    stdshade(t,squeeze(nanmean(allmeas1{c}.naddersp.diff(:,:,2,:),1)-nanmean(allmeas1{c}.naddersp.diff(:,:,1,:))),...
        'r',0.15,1,'std');
    hold on
    stdshade(t,squeeze(nanmean(allmeas2{c}.naddersp.diff(:,:,2,:),1)-nanmean(allmeas2{c}.naddersp.diff(:,:,1,:))),...
        'b',0.15,1,'std');
    Plot_sigmask(p(1,c,1).axis,compstats.ptdiff{c}.prob < 0.05,'cmapline','LineWidth',5);
    xlabel('Time (s)')
    ylabelunits(settings)
    
    FixAxes(gca,14)
    title(settings.tfparams.fbandnames{c})
    
    plotindx = linspace(0,max(settings.aucindex),5);
    plotindx(1) = [];
    for cc = 1:4
        p(1,c,cc+1).select()
        %axes(p(2,c,cc+1).axis)
        plotdata = squeeze(nanmean(allmeas1{c}.naddersp.diff(:,plotindx(cc),2,:)-allmeas1{c}.naddersp.diff(:,plotindx(cc),1,:)-...
            allmeas2{c}.naddersp.diff(:,plotindx(cc),2,:)+allmeas2{c}.naddersp.diff(:,plotindx(cc),1,:),4));
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                0.*compstats.ptdiff{c}.mask(:,plotindx(cc)),0.*compstats.ptdiff{c}.mask(:,plotindx(cc)));
        else
            cluster_topoplot(plotdata,settings.layout,...
                ~(0.*compstats.ptdiff{c}.mask(:,plotindx(cc))),0.*compstats.ptdiff{c}.mask(:,plotindx(cc)));
        end
        colormap(lkcmap2)
        if cc == 4
            cbars1(c) = colorbar('EastOutside');
        end
        ax(cc) = p(1,c,cc+1).axis;
        title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
        Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
    
    % TTV-based differences
    
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    p(2,c).pack();
    for cc = 1:4
        p(2,c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(2,c,1).select();
    stdshade(t,squeeze(nanmean(allmeas1{c}.ttversp.real,1)),...
        'r',0.15,1,'std');
    hold on
    stdshade(t,squeeze(nanmean(allmeas2{c}.ttversp.real,1)),...
        'b',0.15,1,'std');
    Plot_sigmask(p(2,c,1).axis,compstats.ttvdiff{c}.prob < 0.05,'cmapline','LineWidth',5);
    
    FixAxes(gca,14)
    xlabel('Time (s)')
    ylabelunits(settings)
    
    plotindx = linspace(0,max(settings.aucindex),5);
    plotindx(1) = [];
    for cc = 1:4
        p(2,c,cc+1).select()
        %axes(p(2,c,cc+1).axis)
        plotdata = squeeze(nanmean(allmeas1{c}.ttversp.real(:,plotindx(cc),:)-allmeas2{c}.ttversp.real(:,plotindx(cc),:),3));
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                0.*compstats.ttvdiff{c}.mask(:,plotindx(cc)),0.*compstats.ttvdiff{c}.mask(:,plotindx(cc)));
        else
            cluster_topoplot(plotdata,settings.layout,...
                ~(0.*compstats.ttvdiff{c}.mask(:,plotindx(cc))),0.*compstats.ttvdiff{c}.mask(:,plotindx(cc)));
        end
        colormap(lkcmap2)
        if cc == 4
            cbars2(c) = colorbar('EastOutside');
        end
        ax(cc) = p(2,c,cc+1).axis;
        title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
        Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
end

%fix margins

p.margintop = 8;
p.marginleft = 18;
p.marginright = 8;
p(2).de.marginleft = 18;
p(1).de.marginleft = 18;

for c = 1:settings.nfreqs
    ax(c) = p(1,c,1).axis;
    cbars1(c).Position = [ax(c).Position(1)+ax(c).Position(3) ax(c).Position(2) cbars1(c).Position(3) 0.15*ax(c).Position(4)];
    ax(c) = p(2,c,1).axis;
    cbars2(c).Position = [ax(c).Position(1)+ax(c).Position(3) ax(c).Position(2) cbars2(c).Position(3) 0.15*ax(c).Position(4)];
end

AddFigureLabel(p(1,1,1).axis,'A')
AddFigureLabel(p(2,1,1).axis,'B')

set(gcf,'Color','w')

savefig(gcf,'Fig_compare.fig')
export_fig('Fig_compare.png','-m4')
