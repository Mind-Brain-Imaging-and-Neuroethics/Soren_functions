function [stats_mixd,stats_oscifrac,dev_mixd,dev_oscifrac] = Plot_figure_oscifrac_v3(spec,fbands,frange_ple,chanlocs,statsinfo,plotstats,opts)
% Plots a standardized figure for the OsciFrac paper
%
% Input arguments:
%    spec: a structure array. Each field should be titled according to the
%       condition, and should contain an array of IRASA spectra for each
%       subject
%    fbands: an 5x2 array containing the definitions of the frequency bands
%       of interest (the function assumes you're using 5 bands (delta, 
%       theta, alpha, beta, and gamma), but allows you to define the
%       frequency ranges of those bands to deal with different highpasses
%       and lowpasses between datasets)
%    frange_ple: a 1x2 array containing the range over which the PLE should
%       be calculated (estimate this graphically based on when the fractal
%       power "drops off")
%    chanlocs: an EEG chanlocs structure corresponding to the data
%    statfun: the fieldtrip statfun to use for statistics for the data.
%
%    Optional inputs:
%       plotstats: the cluster stats calculated by the function. If this is
%       input, it saves the function from having to calculate them again.



fields = fieldnames(spec);

frange = intersect(find(spec.(fields{1})(1).freq(:,1) > frange_ple(1)),find(spec.(fields{1})(1).freq(:,1) < frange_ple(1,2)));

data = []; data.trial{1} = rand(124,1000); data.time{1} = [1:1000]; data.fsample = 1;
EEG = ft2eeglab(data);

EEG.chanlocs = chanlocs;
EEG.data = rand(length(EEG.chanlocs),1000);
EEG.setname = 'tmpset';
EEG = eeg_checkset(EEG);

data = eeglab2fieldtrip(EEG,'preprocessing','none');
datasetinfo.elec = data.elec;
datasetinfo.label = data.label;


for i = 1:length(fields)
    for c = 1:length(spec.(fields{i}))
        for cc = 1:size(fbands,1)
            %if cc > 1
                relmixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),{'mixd',[fbands(1,1),fbands(end,end)]},'trapz');
            absmixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),0,'mean');
                %else
            %    mixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),0,'mean');
            %end
        end
        
        
        for cc = 1:size(fbands,1)
            %if cc > 1
                relosci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),{'osci',[fbands(1,1),fbands(end,end)]},'trapz');
                absosci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),0,'mean');
                %else
            %    osci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),0,'mean');
            %end
        end
        
        for cc = 1:size(fbands,1)
            %if cc > 1
                %relfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:),{'frac',[fbands(1,1) fbands(end,end)]},'trapz');
                            absfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:),0,'mean');
                %else
            %    frac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:),0,'mean');
            %end
        end
        
        
        %
        %         for cc = 1:6
        %             oscifrac.(fields{i})(c,:,cc) = OsciFrac_EEG_wrapper(spec.(fields{i})(c),fbands(cc,:));
        %         end
        
        tmp = amri_sig_plawfit(spec.(fields{i})(c),frange_ple);
        ple.(fields{i})(c,:) = tmp.Beta;
        
        bb.(fields{i})(c,:) = 1./OsciFrac_EEG_wrapper(spec.(fields{i})(c),[fbands(1,1) fbands(end,end)]);
        %bb.(fields{i})(c,:) = tmp.Cons;
%         spec.(fields{i})(c).linearized_frac = exp(vert(-tmp.Beta)*horz(log(spec.(fields{i})(c).freq))+vert(tmp.Cons));
%         spec.(fields{i})(c).linearized_frac = spec.(fields{i})(c).linearized_frac';
%         spec.(fields{i})(c).resid_frac = spec.(fields{i})(c).frac - spec.(fields{i})(c).linearized_frac;
%         
%         for cc = 1:6
%             if cc > 1
%                 linfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'linearized_frac',fbands(cc,:),0,'mean');
%             else
%                 linfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'linearized_frac',fbands(cc,:),0,'mean');
%             end
%         end
%         
%         for cc = 1:6
%             if cc > 1
%                 resfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'resid_frac',fbands(cc,:),0,'mean');
%             else
%                 resfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'resid_frac',fbands(cc,:),0,'mean');
%             end
%         end
    end
end

for c = 1:length(fields)
    spec.(fields{c}) = mergestructs(spec.(fields{c}));
    spec.(fields{c}) = structfun(@(data)nanmean(data,3),spec.(fields{c}),'UniformOutput',false);
end

fbandnames = {'Delta','Theta','Alpha','Beta','Gamma'};


data = [];
for i = 1:length(statsinfo.statfields)
    tmp = [squeeze(mean(relmixd.(statsinfo.statfields{i}),2)) squeeze(mean(relosci.(statsinfo.statfields{i}),2)) ...
        mean(ple.(statsinfo.statfields{i}),2) mean(bb.(statsinfo.statfields{i}),2) ...
        ones(size(mean(ple.(statsinfo.statfields{i}),2)))*i];
    if strcmpi(statsinfo.paired,'yes') %adding mixed effects
        tmp = cat(2,tmp,1:size(tmp,1));
    end
    data = cat(1,data,tmp);
end

if strcmpi(statsinfo.paired,'yes')
   varnames = [cellcat('Mixd_',fbandnames,'',0) ...
    cellcat('Osci_',fbandnames,'',0) {'PLE' 'Frac_BB' 'Condition' 'Subject'}];
else
   varnames = [cellcat('Mixd_',fbandnames,'',0) ...
    cellcat('Osci_',fbandnames,'',0) {'PLE' 'Frac_BB' 'Condition'}];
end

data = array2table(data,'VariableNames',varnames);

% condn = cellstr(num2str(data.Condition));
% for i = 1:length(statsinfo.statfields)
%     condn = replace(condn,num2str(i),statsinfo.statfields{i});
% end

%data.Condition = condn;

mixnames = cellcat('+',cellcat('Mixd_',fbandnames,'',0),'',1);
mixnames = [mixnames{:}];
mixnames(end) = [];

oscinames = cellcat('+',cellcat('Osci_',fbandnames,'',0),'',1);
oscinames = [oscinames{:}];

frmla1 = ['Condition ~ ' mixnames];
frmla2 = ['Condition ~ ' oscinames 'PLE+Frac_BB'];

if strcmpi(statsinfo.paired, 'yes')
   frmla1 = [frmla1 '+(1|Subject)-1'];
   frmla2 = [frmla2 '+(1|Subject)-1'];
end

[~,dev_mixd,stats_mixd] = mnrfit(zscore(data{:,1:size(fbands,1)},[],1),data{:,end},'model','nominal');

[~,dev_oscifrac,stats_oscifrac] = mnrfit(zscore(data{:,(1+size(fbands,1)):(7+size(fbands,1))},[],1),data{:,end},'model','nominal');


p = panel('no-manage-font');

p.pack('v',{40 60})

p(1).select()
for c = 1:5
    for cc = 1:length(fields)
        plotstack(c,cc,1) = nanmean(nanmean(absfrac.(fields{cc})(:,:,c),2),1);
        plotstack(c,cc,2) = nanmean(nanmean(absosci.(fields{cc})(:,:,c),2),1);
        %plotstack(c,cc,2) = nanmean(nanmean(linfrac.(fields{cc})(:,:,c),2),1);
        %plotstack(c,cc,3) = nanmean(nanmean(resfrac.(fields{cc})(:,:,c),2),1);
    end
end

l = lines;
l = l(1:3,:);
for c = 1:3
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

[h,barpos] = plotBarStackGroups(plotstack,fbandnames,col);
hold on;
for cc = 1:length(fields)
    if size(absmixd.(fields{cc}),1) > 1
        er = errorbar(barpos(cc,:),squeeze(nanmean(nanmean(absmixd.(fields{cc}),1),2)),squeeze(nanstd(nanmean(absmixd.(fields{cc}),2),[],1)));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
end
legend(l)
FixAxes(gca,14)

p(2).pack('h',{35 45 20})
pris = prism;
pris = pris(1:6,:);
plot_mixnames = {'Mixed Delta','Mixed Theta','Mixed Alpha','Mixed Beta','Mixed Gamma'};
p(2,1).select()
for c = 1:size(stats_mixd.beta,2)
    scatter((1:size(stats_mixd.beta(2:end,c),1))+c*0.15-0.15,stats_mixd.beta(2:end,c),36,pris(c,:));
    hold on
    er = errorbar((1:size(stats_mixd.beta(2:end,c),1))+c*0.15-0.15,stats_mixd.beta(2:end,c),stats_mixd.se(2:end,c),'LineStyle','none','Color',pris(c,:),'LineWidth',2);
end
xax = get(gca,'XLim');
line(xax,[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5])
set(gca,'XLim',xax+[-0.1 0.1],'XTick',(1:size(stats_mixd.beta(2:end,:),1))+(size(stats_mixd.beta,2)*0.15-0.15)/2,'XTickLabel',plot_mixnames)
ylabel('Coefficient')
FixAxes(gca,14)

plot_oscifracnames = {'Osci Delta','Osci Theta','Osci Alpha','Osci Beta','Osci Gamma','Frac PLE','Frac broadband'};
p(2,2).select()
for c = 1:size(stats_mixd.beta,2)
    scatter((1:size(stats_oscifrac.beta(2:end,c),1))+c*0.15-0.15,stats_oscifrac.beta(2:end,c),36,pris(c,:));
    hold on
    er = errorbar((1:size(stats_oscifrac.beta(2:end,c),1))+c*0.15-0.15,stats_oscifrac.beta(2:end,c),stats_oscifrac.se(2:end,c),'LineStyle','none','Color',pris(c,:),'LineWidth',2);
end
xax = get(gca,'XLim');
line(xax,[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5])
set(gca,'XLim',xax+[-0.1 0.1],'XTick',(1:length(stats_oscifrac.beta(2:end)))+(size(stats_oscifrac.beta,2)*0.15-0.15)/2,'XTickLabel',plot_oscifracnames)
ylabel('Coefficient')
FixAxes(gca,14)

p(2,3).select()
se_dev = ones(2,1); % figure out how to get errors or CIs for deviance, AIC, etc
er = errorbar([dev_mixd dev_oscifrac],se_dev,'LineStyle','none','LineWidth',2,'Color',[0 0 0]);
set(gca,'XTickLabel',{'Mixed power','IRASA decomposition'},'XLim',get(gca,'XLim')+[-0.1 0.1])


