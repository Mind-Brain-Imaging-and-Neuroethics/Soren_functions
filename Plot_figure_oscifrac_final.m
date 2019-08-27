function [stats_mixd,stats_oscifrac,cmpare,newcmpare] = Plot_figure_oscifrac_bayes(spec,fbands,opts)
% Plots a standardized figure for the OsciFrac paper
%
% Input arguments:
%    spec: a structure array. Each field should be titled according to the
%       condition, and should contain an array of IRASA spectra for each
%       subject
%    fbands: an nx2 array containing the definitions of the frequency bands
%       of interest
%    opts: a structure with the following fields:
%       paired: 'yes' or 'no' for paired statistics. Currently not used
%          (default = 'no')
%       statfields: fields to do statistics on (default = all the fields in
%          the specs structure). The order of this field also specifies the
%          reference category - the last field is taken as the reference
%       fbandnames: names of the frequency bands you're using (default =
%          {'Delta','Theta','Alpha','Beta','Gamma'})
%       frange_ple: a 1x2 array containing the range over which the PLE
%          should be calculated (estimate this graphically based on when
%          the fractal power "drops off"; default = [lowest frequency in
%          fbands highest frequency in fbands])
%       stats_mixd: can input the stats from a previous call to this
%          function in order to do the plot faster
%       stats_oscifrac: can input the stats from a previous call to this
%          function in order to do the plot faster
%       cmpare: can input the stats from a previous call to this
%          function in order to do the plot faster

%% Setting up options
fields = fieldnames(spec);

if nargin < 3
    opts = struct;
end

opts = setdefault(opts,'paired','no');
opts = setdefault(opts,'fbandnames',{'Delta','Theta','Alpha','Beta','Gamma'});
opts = setdefault(opts,'statfields',fields);
opts = setdefault(opts,'frange_ple',[fbands(1,1) fbands(end,end)]);

frange = intersect(find(spec.(fields{1})(1).freq(:,1) > opts.frange_ple(1)),find(spec.(fields{1})(1).freq(:,1) < opts.frange_ple(1,2)));

%% Calculating power
for i = 1:length(fields)
    for c = 1:length(spec.(fields{i}))
        for cc = 1:size(fbands,1)
            relmixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),{'mixd',[fbands(1,1),fbands(end,end)]},'trapz');
            absmixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),0,'mean');
        end
        
        for cc = 1:size(fbands,1)
            relosci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),{'osci',[fbands(1,1),fbands(end,end)]},'trapz');
            absosci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),0,'mean');
        end
        
        for cc = 1:size(fbands,1)
            absfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:),0,'mean');
        end
        
        tmp = amri_sig_plawfit(spec.(fields{i})(c),opts.frange_ple);
        ple.(fields{i})(c,:) = tmp.Beta;
        
        bb.(fields{i})(c,:) = 1./OsciFrac_EEG_wrapper(spec.(fields{i})(c),[fbands(1,1) fbands(end,end)]);
    end
end

% for c = 1:length(fields)
%     spec.(fields{c}) = mergestructs(spec.(fields{c}));
%     spec.(fields{c}) = structfun(@(data)nanmean(data,3),spec.(fields{c}),'UniformOutput',false);
% end

%% Logistic regression

if ~isfield(opts,'stats_mixd') && ~isfield(opts,'stats_oscifrac') && ~isfield(opts,'cmpare')
    
    data = [];
    for i = 1:length(opts.statfields)
        tmp = [squeeze(mean(absmixd.(opts.statfields{i}),2)) squeeze(mean(absosci.(opts.statfields{i}),2)) ...
            mean(ple.(opts.statfields{i}),2) mean(bb.(opts.statfields{i}),2) ...
            ones(size(mean(ple.(opts.statfields{i}),2)))*i];
        data = cat(1,data,tmp);
    end
    
    if strcmpi(opts.paired,'yes') %adding subjects for random effects
        data = cat(2,data,repmat([1:size(tmp,1)]',length(fields),1));
    end
    
    if strcmpi(opts.paired,'yes')
        varnames = [cellcat('Mixd_',opts.fbandnames,'',0) ...
            cellcat('Osci_',opts.fbandnames,'',0) {'PLE' 'Frac_BB' 'Condition' 'Subject'}];
    else
        varnames = [cellcat('Mixd_',opts.fbandnames,'',0) ...
            cellcat('Osci_',opts.fbandnames,'',0) {'PLE' 'Frac_BB' 'Condition'}];
    end
    
    data = array2table(data,'VariableNames',varnames);
    
    mixnames = cellcat('+',cellcat('Mixd_',opts.fbandnames,'',0),'',1);
    mixnames = [mixnames{:}];
    mixnames(end) = [];
    
    oscinames = cellcat('+',cellcat('Osci_',opts.fbandnames,'',0),'',1);
    oscinames = [oscinames{:}];
    
    frmla1 = ['Condition ~ ' mixnames];
    frmla2 = ['Condition ~ ' oscinames 'PLE+Frac_BB'];
    
    if strcmpi(opts.paired, 'yes')
        frmla1 = [frmla1 '+(1|Subject)-1'];
        frmla2 = [frmla2 '+(1|Subject)-1'];
    end
    
    if strcmpi(opts.paired,'yes')
        data{:,1:end-2} = zscore(data{:,1:end-2},[],1);
    else
        data{:,1:end-1} = zscore(data{:,1:end-1},[],1);
    end
    
    if strcmpi(opts.paired,'yes')
        for c = 1:size(fbands,1)
            t = rm_anova2(cat(1,data.(['Mixd_' opts.fbandnames{c}]),data.(['Osci_' opts.fbandnames{c}])),...
                cat(1,data.Subject,data.Subject),cat(1,data.Condition,data.Condition),...
                Make_designVect([height(data) height(data)])',{'Condition','Oscifrac'});
            p_dif_irasa(c) = t{4,6};
        end
    else
        for c = 1:size(fbands,1)
            [~,t] = anovan(cat(1,data.(['Mixd_' opts.fbandnames{c}]),data.(['Osci_' opts.fbandnames{c}])),...
                {Make_designVect([height(data) height(data)])',cat(1,data.Condition,data.Condition)},'model','interaction');
            p_dif_irasa(c) = t{4,7};
            %get p value for interaction term here
        end
    end

    %[stats_mixd,stats_oscifrac,cmpare] = mxmnr_compare(data,frmla1,frmla2);
    
    
    
    
%     
%     %stats_mixd.fixed = mergestructs(stats_mixd.fixed);
%     %stats_oscifrac.fixed = mergestructs(stats_oscifrac.fixed);
%     
%     stats_mixd.beta = extractfield(stats_mixd.fixed,'Estimate');
%     stats_mixd.ci = [extractfield(stats_mixd.fixed,'l_95_CI')' extractfield(stats_mixd.fixed,'u_95_CI')'];
%     stats_oscifrac.beta = extractfield(stats_oscifrac.fixed,'Estimate');
%     stats_oscifrac.ci = [extractfield(stats_oscifrac.fixed,'l_95_CI')' extractfield(stats_oscifrac.fixed,'u_95_CI')'];
%     
%     if ~strcmpi(opts.paired,'yes')
%         stats_mixd.beta(1) = [];
%         stats_oscifrac.beta(1) = [];
%         stats_mixd.ci(1,:) = [];
%         stats_oscifrac.ci(1,:) = [];
%     end
%     
%     stats_mixd.beta = reshape(stats_mixd.beta,length(opts.fbandnames),[]);
%     stats_oscifrac.beta = reshape(stats_oscifrac.beta,length(opts.fbandnames)+2,[]);
%     stats_mixd.ci = cat(2,reshape(stats_mixd.ci(:,1),length(opts.fbandnames),1,[]),...
%         reshape(stats_mixd.ci(:,2),length(opts.fbandnames),1,[]));
%     stats_oscifrac.ci = cat(2,reshape(stats_oscifrac.ci(:,1),length(opts.fbandnames)+2,1,[]),...
%         reshape(stats_oscifrac.ci(:,2),length(opts.fbandnames)+2,1,[]));
%     
%     
%     
%     newcmpare.elpd_mixd = stats_mixd.kfold.estimates(1).Estimate;
%     newcmpare.elpd_oscifrac = stats_oscifrac.kfold.estimates(1).Estimate;
%     newcmpare.elpd_mixd_se = stats_mixd.kfold.estimates(1).SE;
%     newcmpare.elpd_oscifrac_se = stats_oscifrac.kfold.estimates(1).SE;
%     newcmpare.se_diff = cmpare(2).se_diff;
else
    stats_mixd = opts.stats_mixd;
    stats_oscifrac = opts.stats_oscifrac;
    newcmpare = opts.newcmpare;
    cmpare = opts.cmpare;
end

%[~,~,stats_mixd] = mnrfit(zscore(data{:,1:size(fbands,1)},[],1),data.Condition,'model','nominal');

%[~,~,stats_oscifrac] = mnrfit(zscore(data{:,(1+size(fbands,1)):(7+size(fbands,1))},[],1),data.Condition,'model','nominal');

%% Plotting the figure
p = panel('no-manage-font');

p.pack('v',{40 60})

p(1).select()
for c = 1:size(fbands,1)
    for cc = 1:length(fields)
        plotstack(c,cc,1) = nanmean(nanmean(absfrac.(fields{cc})(:,:,c),2),1);
        plotstack(c,cc,2) = nanmean(nanmean(absosci.(fields{cc})(:,:,c),2),1);
        mixsum(c,cc) = nanmean(nanmean(absmixd.(fields{cc})(:,:,c),2),1);
    end
end

l = lines;
l = l(1:length(fields),:);
for c = 1:length(fields)
    palel(c,:) = palecol(l(c,:));
end
for c = 1:size(fbands,1)
    col{c} = cat(3,l,palel);
    col{c} = permute(col{c},[3 1 2]);
end

l = [];
for c = 1:length(fields)
    l = [l {[fields{c} ' Fractal']} {[fields{c} ' Oscillatory']}];
end

[h,barpos] = plotBarStackGroups(plotstack,opts.fbandnames,col);
hold on;
for cc = 1:length(fields)
    if size(absmixd.(fields{cc}),1) > 1
        er = errorbar(barpos(cc,:),mixsum(:,cc),squeeze(nanstd(nanmean(absmixd.(fields{cc}),2),[],1)));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
end
legend(l)
FixAxes(gca,14)
ylabel('Power')
set(gca,'YScale','log')

p(2).pack('h',{35 45 20})
%pris = prism;
%pris = pris(5:10,:);
pris = lines;
pris = pris(2:7,:);
plot_mixnames = cellcat('Mixed ',opts.fbandnames,'',0);
p(2,1).select()
for c = 1:size(stats_mixd.beta,2)
    scatter((1:size(stats_mixd.beta(:,c),1))+c*0.15-0.15,stats_mixd.beta(:,c),36,pris(c,:),'HandleVisibility','off');
    hold on
    er = errorbar((1:size(stats_mixd.beta(:,c),1))+c*0.15-0.15,stats_mixd.beta(:,c),...
        stats_mixd.ci(:,1,c)-stats_mixd.beta(:,c),abs(stats_mixd.ci(:,2,c)-stats_mixd.beta(:,c)),...
        'LineStyle','none','Color',pris(c,:),'LineWidth',2);
end
xax = get(gca,'XLim');
line(xax,[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5])
set(gca,'XLim',xax+[-0.1 0.1],'XTick',...
    (1:size(stats_mixd.beta(:,:),1))+(size(stats_mixd.beta,2)*0.15-0.15)/2,'XTickLabel',plot_mixnames)
ylabel('Coefficient')
legend(cellcat('Likelihood ',cellcat(opts.statfields{1},opts.statfields(2:end),'-',0),'',0))
FixAxes(gca,14)
fix_xticklabels(gca,0.1,{'FontSize',12})

plot_oscifracnames = [cellcat('Osci ',opts.fbandnames,'',0) {'Frac PLE','Frac broadband'}];
p(2,2).select()
for c = 1:size(stats_mixd.beta,2)
    scatter((1:size(stats_oscifrac.beta(:,c),1))+c*0.15-0.15,stats_oscifrac.beta(:,c),36,pris(c,:),'HandleVisibility','off');
    hold on
    er = errorbar((1:size(stats_oscifrac.beta(:,c),1))+c*0.15-0.15,stats_oscifrac.beta(:,c),...
        stats_oscifrac.ci(:,1,c)-stats_oscifrac.beta(:,c),...
        abs(stats_oscifrac.ci(:,2,c)-stats_oscifrac.beta(:,c)),...
        'LineStyle','none','Color',pris(c,:),'LineWidth',2);
end
xax = get(gca,'XLim');
line(xax,[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5])
set(gca,'XLim',xax+[-0.1 0.1],'XTick',...
    (1:length(stats_oscifrac.beta(2:end)))+(size(stats_oscifrac.beta,2)*0.15-0.15)/2,'XTickLabel',plot_oscifracnames)
ylabel('Coefficient')
legend(cellcat('Likelihood ',cellcat(opts.statfields{1},opts.statfields(2:end),'-',0),'',0))
FixAxes(gca,14)
fix_xticklabels(gca,0.1,{'FontSize',12})

p(2,3).select()
scatter([1,2],[newcmpare.elpd_mixd newcmpare.elpd_oscifrac],72,[0 0 0],'x','LineWidth',2);
hold on
errorbar([1 2],[newcmpare.elpd_mixd newcmpare.elpd_oscifrac],4*[newcmpare(1).elpd_mixd_se newcmpare(1).elpd_oscifrac_se],'LineStyle','none','Color',[0 0 0],'LineWidth',2)
set(gca,'XLim',get(gca,'XLim')+[-0.1 0.1])
set(gca,'XTick',[1 2],'XTickLabel',{'Mixed power','IRASA decomposition'})
ylabel('Expected Log Predictive Density')
FixAxes(gca,14)
fix_xticklabels(gca,0.1,{'FontSize',12})

p.marginleft = 24;
p.marginbottom = 18;

