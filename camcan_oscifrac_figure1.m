function camcan_oscifrac_figure1(spec,fbands,opts)

fields = fieldnames(spec);

if nargin < 3
   opts = struct; 
end

opts = setdefault(opts,'fbandnames',{'Delta','Theta','Alpha','Beta','Gamma'});
opts = setdefault(opts,'statfields',fields);
opts = setdefault(opts,'frange_ple',[fbands(1,1) fbands(end,end)]);

frange = intersect(find(spec(1).freq(:,1) > opts.frange_ple(1)),find(spec(1).freq(:,1) < opts.frange_ple(1,2)));

%% Calculating power
for c = 1:length(spec)
    for cc = 1:size(fbands,1)
        relmixd(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'mixd',fbands(cc,:),{'mixd',[fbands(1,1),fbands(end,end)]},'trapz');
        absmixd(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'mixd',fbands(cc,:),0,'mean');
    end
    
    for cc = 1:size(fbands,1)
        relosci(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'osci',fbands(cc,:),{'osci',[fbands(1,1),fbands(end,end)]},'trapz');
        absosci(c,:,cc) = IRASAPower_EEG_wrapper(spec(c),'osci',fbands(cc,:),0,'mean');
    end
    
    %         for cc = 1:size(fbands,1)
    %             absfrac.(c,:,cc) = IRASAPower_EEG_wrapper(spec.(c),'frac',fbands(cc,:),0,'mean');
    %         end
    
    for cc = 1:size(fbands,1)
       oscifrac(c,:,cc) = OsciFrac_EEG_wrapper(spec(c),fbands(cc,:));
    end
    
    tmp = amri_sig_plawfit(spec.(c),opts.frange_ple);
    ple(c,:) = tmp.Beta;
    
    bb(c,:) = 1./OsciFrac_EEG_wrapper(spec.(c),[fbands(1,1) fbands(end,end)]);
end

meanspec.(fields{c}) = mergestructs(spec.(fields{c}));
meanspec.(fields{c}) = structfun(@(data)nanmean(data,3),meanspec.(fields{c}),'UniformOutput',false);

%% Regression stuff

for i = 1:size(fbands,1)
    tbl = array2table([nanmean(relmixd(:,:,i),2) nanmean(relosci(:,:,i),2) nanmean(ple,2) nanmean(bb,2)],...
        'VariableNames',{['Mixed_' opts.fbandnames{i}] ['Osci_' opts.fbandnames{i}] 'PLE' 'Frac_BB'});
   mdl{i} = fitlm(tbl,['Mixed_' opts.fbandnames{i} ' ~ Osci_' opts.fbandnames{i} ' + PLE + Frac_BB + PLE*Frac_BB']);
end


%% Plotting the figure

figure

p = panel('no-manage-font');

p.pack('h',{30 70})

p(1,1).select()
hold on
loglog(meanspec.freq(frange,1),nanmean(meanspec.mixd(frange,:),2),'LineWidth',1.5)

legend(fields)
FixAxes(gca,12)
title('Mixed power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XScale','log','YScale','log','XLim',opts.frange_ple)

p(1,2,1).select()
hold on
plot(meanspec.freq(frange,1),nanmean(meanspec.osci(frange,:),2),'LineWidth',1.5)

FixAxes(gca,12)
legend(fields)
title('Oscillatory power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XLim',opts.frange_ple)

p(1,2,2).select()
hold on
loglog(meanspec.freq(frange,1),nanmean(meanspec.frac(frange,:),2),'LineWidth',1.5)
legend(fields)
FixAxes(gca,12)
title('Fractal power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XScale','log','YScale','log','XLim',opts.frange_ple)

p(2).pack(2,size(fbands,1))

for i = 1:size(fbands,1)
   p(2,1,i).select()
   ft_topoplot_vec('neuromag306mag.lay',mean(oscifrac,1),opts.datasetinfo.label)
   
   p(2,2,i).select()
   scatter(1:3,mdl{i}.Coefficients.Estimate(2:end))
   er = errorbar(1:3,mdl{i}.Coefficients.Estimate(2:end),mdl{i}.Coefficients.SE*1.96,...
       'LineStyle','none','LineWidth',2);
   xl = get(gca,'XLim');
   line(xl,[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5])
   set(gca,'XLim',xl,'XTickLabel',{['Osci_' opts.fbandnames{i}],'PLE','Frac_BB'})
   fix_xticklabels(gca)
end





