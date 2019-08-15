function [data] = camcan_oscifrac_figure1(spec,fbands,opts)

fields = fieldnames(spec);

if nargin < 3
    opts = struct;
end

opts = setdefault(opts,'fbandnames',{'Delta','Theta','Alpha','Beta','Gamma'});
opts = setdefault(opts,'statfields',fields);
opts = setdefault(opts,'frange_ple',[fbands(1,1) fbands(end,end)]);



frange = intersect(find(spec(1).freq(:,1) > opts.frange_ple(1)),find(spec(1).freq(:,1) < opts.frange_ple(1,2)));

%% Calculating power

if ~isfield(opts,'data')
    for c = 1:length(spec)
        for cc = 1:size(fbands,1)
            relmixd(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'mixd',fbands(cc,:),{'mixd',[fbands(1,1),fbands(end,end)]},'trapz',1);
            absmixd(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'mixd',fbands(cc,:),0,'mean',1);
        end
        
        for cc = 1:size(fbands,1)
            relosci(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'osci',fbands(cc,:),{'osci',[fbands(1,1),fbands(end,end)]},'trapz',1);
            absosci(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'osci',fbands(cc,:),0,'mean',1);
        end
        
        for cc = 1:size(fbands,1)
            relfrac(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'frac',fbands(cc,:),{'frac',[fbands(1,1),fbands(end,end)]},'trapz',1);
            absfrac(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'frac',fbands(cc,:),0,'mean');
        end
        
        for cc = 1:size(fbands,1)
            oscifrac(c,:,cc) = OsciFrac_EEG_wrapper(spec(c),fbands(cc,:));
        end
        
        tmp = amri_sig_plawfit(spec(c),opts.frange_ple);
        ple(c,:) = tmp.Beta;
        
        bb(c,:) = 1./OsciFrac_EEG_wrapper(spec(c),[fbands(1,1) fbands(end,end)]);
    end
else
    datafields = fieldnames(opts.data);
    for c = 1:length(datafields)
        eval([datafields{c} '= opts.data.(datafields{c});']);
    end
end

meanspec = mergestructs(spec);
meanspec = structfun(@(data)nanmean(data,3),meanspec,'UniformOutput',false);

%% Regression stuff
for i = 1:size(fbands,1)
    tbl = array2table([nanmean(absmixd(:,:,i),2) nanmean(absosci(:,:,i),2) nanmean(ple,2) nanmean(bb,2)],...
        'VariableNames',{['Mixed_' opts.fbandnames{i}] ['Osci_' opts.fbandnames{i}] 'PLE' 'Frac_BB'});
    tbl{:,:} = zscore(tbl{:,:},[],1);
    mdl{i} = fitlm(tbl,['Mixed_' opts.fbandnames{i} ' ~ Osci_' opts.fbandnames{i} ' + PLE + Frac_BB']);
    
    partr(i,1) = semipartcorr_resid(tbl{:,2},tbl{:,1},tbl{:,3:4});
    partr(i,2) = semipartcorr_resid(tbl{:,3},tbl{:,1},tbl{:,[2 4]});
    partr(i,3) = semipartcorr_resid(tbl{:,4},tbl{:,1},tbl{:,[2 3]});
    
    for c = 1:3
       partr_ci(i,c,:) = corrci(partr(i,c),height(tbl)); 
    end
end

%partr = partr.^2;
%partr_ci = partr_ci.^2;

%% Plotting the figure

figure

p = panel('no-manage-font');

p.pack('v',{50 50})

p(1).pack('v',{50 50})
p(1,2).pack('h',{50 50})

p(1,1).select()
hold on
loglog(meanspec.freq(frange,1),nanmean(meanspec.mixd(frange,:),2),'LineWidth',2)

FixAxes(gca,14)
title('Mixed power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XScale','log','YScale','log','XLim',opts.frange_ple)

p(1,2,1).select()
hold on
plot(meanspec.freq(frange,1),nanmean(meanspec.osci(frange,:),2),'LineWidth',2)

FixAxes(gca,14)
title('Oscillatory power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XLim',opts.frange_ple)

p(1,2,2).select()
hold on
loglog(meanspec.freq(frange,1),nanmean(meanspec.frac(frange,:),2),'LineWidth',2)
FixAxes(gca,14)
title('Fractal power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XScale','log','YScale','log','XLim',opts.frange_ple)


p.marginleft = 22;
p(1).marginbottom = 32;
p(1,1).marginbottom = 24;
p.margintop = 8;

p(2).pack(2,size(fbands,1))

count = 1;
for i = 1:size(fbands,1)
    p(2,1,i).select()
    ft_topoplot_vec('neuromag306mag.lay',mean(absosci(:,:,i),1),opts.datasetinfo.label)
    title(['Oscillatory power' newline opts.fbandnames{i}],'FontSize',14)
    colorbar
    %ax(count) = gca;
    %count = count+1;
    
    p(2,2,i).select()
    ft_topoplot_vec('neuromag306mag.lay',mean(absfrac(:,:,i),1),opts.datasetinfo.label)
    title(['Fractal power' newline opts.fbandnames{i}],'FontSize',14)
    %if i == size(fbands,1)
        colorbar
    %end
    %ax(count) = gca;
    %count = count+1;
end
%Normalize_Clim(ax,0)
%clear ax


figure

p = panel('no-manage-font')

p.pack('v',{1/4 1/4 1/4 1/4})


p.marginright = 10;

p(1).select()

for c = 1:size(fbands,1)
    plotstack(c,1,1) = nanmean(nanmean(absfrac(:,:,c),2),1);
    plotstack(c,1,2) = nanmean(nanmean(absosci(:,:,c),2),1);
end

l = lines;
l = l(1,:);
for c = 1:1
    palel(c,:) = palecol(l(c,:));
end
for c = 1:size(fbands,1)
    col{c} = cat(3,l,palel);
    col{c} = permute(col{c},[3 1 2]);
end

[h,barpos] = plotBarStackGroups(plotstack,opts.fbandnames,col);
hold on;
for cc = 1:1
    if size(absmixd,1) > 1
        er = errorbar(barpos(cc,:),squeeze(nanmean(nanmean(absmixd,1),2)),squeeze(nanstd(nanmean(absmixd,2),[],1)));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
end
legend({'Fractal Power','Oscillatory Power'})
FixAxes(gca,14)
ylabel('Power')

p(2).pack('h',repmat({1/size(fbands,1)},1,size(fbands,1)))
for i = 1:size(fbands,1)
    p(2,i).select()
    ft_topoplot_vec('neuromag306mag.lay',mean(oscifrac(:,:,i),1),opts.datasetinfo.label)
    title(['Oscillatory-Fractal ratio' newline opts.fbandnames{i}],'FontSize',14)
    %if i == size(fbands,1)
    colorbar
    %end
    %ax(i) = gca;
end
%Normalize_Clim(ax,0)

p(3).select()

l = lines;
for i = 1:size(fbands,1)
    %scatter((1:3)+0.1*i,mdl{i}.Coefficients.Estimate(2:end),36,l(i,:),'x','LineWidth',2)
    hold on
    er = errorbar((1:3)+0.1*i,mdl{i}.Coefficients.Estimate(2:end),mdl{i}.Coefficients.SE(2:end)*1.96,...
        'LineStyle','none','LineWidth',2,'Color',l(i,:));
end
FixAxes(gca,14)
xl = get(gca,'XLim');
line(xl,[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off')
set(gca,'XLim',xl,'XTick',(1:3)+0.1*(size(fbands,1)/2+1),...
    'XTickLabel',{'Oscillatory Power','Fractal PLE','Fractal Broadband'})
legend(opts.fbandnames)
ylabel('Regression coefficient')


p(4).pack('h',repmat({1/size(fbands,1)},1,size(fbands,1)))
for i = 1:size(fbands,1)
    p(4,i).select()
    hold on
    scatter(1:3,partr(i,:),36,[0 0 0],'o','LineWidth',2)
    %er = errorbar(1:3,partr(i,:),partr(i,:)-partr_ci(i,:,1),partr_ci(i,:,2)-partr(i,:),...
    %    'LineStyle','none','LineWidth',2,'Color',[0 0 0],'HandleVisibility','off');

    xl = get(gca,'XLim');
    line(xl+[-0.1 0.1],[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
    set(gca,'XLim',xl + [-0.1 0.1],'XTickLabel',{['Osci ' opts.fbandnames{i}],'PLE','Frac BB'},'YLim',[-1 1])
    %legend({'Regression Coefficient','Partial R^2'})
    if i == 1
    ylabel('Partial R^2')
    end
    FixAxes(gca,14)
    fix_xticklabels(gca,0.1,{'FontSize',14})
end

p.marginleft = 20;
p(1).marginbottom = 24;

data.absosci = absosci;
data.absmixd = absmixd;
data.absfrac = absfrac;
data.relosci = relosci;
data.relfrac = relfrac;
data.relmixd = relmixd;
data.ple = ple;
data.bb = bb;
data.oscifrac = oscifrac;



