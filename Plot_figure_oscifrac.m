function [p,plotstats] = Plot_figure_oscifrac(spec,fbands,frange_ple,chanlocs,statfun,plotstats,opts)
% Plots a standardized figure for the OsciFrac paper
%
% Input arguments:
%    spec: a structure array. Each field should be titled according to the
%       condition, and should contain an array of IRASA spectra for each
%       subject
%    fbands: an 6x2 array containing the definitions of the frequency bands
%       of interest (the function assumes you're using 6 bands (broadband,
%       delta, theta, alpha, beta, and gamma), but allows you to define the
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
%         for cc = 1:6
%             if cc > 1
%                 mixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),{'mixd',fbands(1,:)},'mean');
%             else
%                 mixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:));
%             end
%         end
        
        for cc = 1:6
            if cc > 1
                osci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),{'osci',fbands(1,:)},'mean');
            else
                osci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:));
            end
        end
        
        for cc = 1:6
            if cc > 1
                frac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:),{'frac',fbands(1,:)},'mean');
            else
                frac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:));
            end
        end
        
%         for cc = 1:6
%             oscifrac.(fields{i})(c,:,cc) = OsciFrac_EEG_wrapper(spec.(fields{i})(c),fbands(cc,:));
%         end
        
        ples.(fields{i})(c,:,cc) = PLE_IRASA_EEG_wrapper(spec.(fields{i})(c),frange_ple); 
    end
end

for c = 1:length(fields)
spec.(fields{c}) = mergestructs(spec.(fields{c}));
spec.(fields{c}) = structfun(@(data)nanmean(data,3),spec.(fields{c}),'UniformOutput',false);
end


p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*3 pos(4)*3],'color','w');

p.de.margin = [5 5 5 5];

%fix margins
p.marginleft = 18;
p.margintop = 8;

p.pack('h',{1/3 2/3})

p(1).pack('v',{1/2 1/2})
p(1,1).marginbottom = 24;


p(1,2).pack('h',{1/2 1/2})

cmap = lines;

p(1,1).select()
hold on
for c = 1:length(fields)
loglog(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).mixd(frange,:),2),'LineWidth',1.5)
end
legend(fields)
FixAxes(gca,12)
title('Mixed power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XScale','log','YScale','log','XLim',frange_ple)

p(1,2,1).select()
hold on
for c = 1:length(fields)
plot(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).osci(frange,:),2),'LineWidth',1.5)
end
FixAxes(gca,12)
legend(fields)
title('Oscillatory power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XLim',frange_ple)

p(1,2,2).select()
hold on
for c = 1:length(fields)
    loglog(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).frac(frange,:),2),'LineWidth',1.5)
end
legend(fields)
FixAxes(gca,12)
title('Fractal power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XScale','log','YScale','log','XLim',frange_ple)


p(2).pack('v',{50 50})
%p(2,1).pack(2,6)
%p(2,2).pack(2,6)
p(2,1).pack(2,3)
p(2,1).de.margin = [10 10 5 5];
p(2,2).pack(2,3)
p(2,2).de.margin = [10 10 5 5];
%p(2,2).pack('h',{1/6 1/6 1/6 1/6 1/6 1/6})

plotstruct = struct;
for c = 1:length(fields)
    plotstruct.(fields{c}) = [];
end

fbandnames = {'Broadband','Delta','Theta','Alpha','Beta','Gamma'};
%ple_fbands = {'2-8 Hz','8-13 Hz','13-50 Hz'};
%ple_fbands = {'2-50 Hz'};

if ~exist('plotstats','var')
plotstats = cell(2,6);
end

if ~strcmpi(statfun,'ft_statfun_friedman') && ~strcmpi(statfun,'ft_statfun_kruskal')
   numstat = 2;
else
    numstat = length(fields);
end

for c = 1:6
    %p(2,1,floor(c/4)+1,c-floor(c/4)*3).select()
    x = floor(c/4)+1;
    y = c-floor(c/4)*3;
    p(2,1,x,y).pack()
   
    
    p(2,1,x,y,1).select()
    if x == 1 && y == 1
       AddFigureLabel(p(2,1,1,1,1).axis,'Osci','middle_left','FontSize',14)
    end
    
    for cc = 1:length(fields)
        plotstruct.(fields{cc}) = squeeze(nanmean(osci.(fields{cc})(:,:,c),2));
    end
    violinplot(plotstruct)
    
    yax = get(gca,'YLim');
    
    set(gca,'YLim',[yax(1) yax(2)+(0.4*(yax(2)-yax(1)))]);
    
    for cc = 1:length(fields)
       p(2,1,x,y).pack({[(cc-1)*1/length(fields) 0.7 1/length(fields) 0.3]})  
    end
    
    title(fbandnames{c})    
    %p(2,1,c,2).pack('h',repmat({1/length(fields)},1,length(fields)))
    for cc = 1:(min(length(fields),numstat))
        statdata{cc} = osci.(fields{cc})(:,:,c)';
    end
    
    opts.nrand = 1000;
    opts.minnbchan = 0;
    %opts.parpool = 24; %remove later
    if isempty(plotstats{1,c})
        plotstats{1,c} = EasyClusterCorrect(statdata,datasetinfo,statfun,opts);
    end
    
    
    for cc = 1:length(fields)
        p(2,1,x,y,cc+1).select()
        cluster_topoplot(nanmean(osci.(fields{cc})(:,:,c),1),chanlocs,plotstats{1,c}.prob,plotstats{1,c}.prob < 0.05,1);
        ax(cc) = gca;
        if cc == length(fields)
           pos = ax(cc).Position;
           cbar = colorbar('EastOutside');
           %set(cbar,'Position',[pos(1)+pos(3)-cbar.Position(3) cbar.Position(2:4)]);
        end
    end
    Normalize_Clim(ax,0)
end

for c = 1:6
    %p(2,1,floor(c/4)+1,c-floor(c/4)*3).select()
    x = floor(c/4)+1;
    y = c-floor(c/4)*3;
    p(2,2,x,y).pack()
    
    p(2,2,x,y,1).select()
        if x == 1 && y == 1
       AddFigureLabel(p(2,2,1,1,1).axis,'Frac','middle_left','FontSize',14)
    end
    
    for cc = 1:length(fields)
        plotstruct.(fields{cc}) = squeeze(nanmean(frac.(fields{cc})(:,:,c),2));
    end
    violinplot(plotstruct)
    
    yax = get(gca,'YLim');
    
    set(gca,'YLim',[yax(1) yax(2)+(0.4*(yax(2)-yax(1)))]);
    
    for cc = 1:length(fields)
       p(2,2,x,y).pack({[(cc-1)*1/length(fields) 0.7 1/length(fields) 0.3]})  
    end
    
    title(fbandnames{c})    
    %p(2,1,c,2).pack('h',repmat({1/length(fields)},1,length(fields)))
    for cc = 1:(min(length(fields),numstat))
        statdata{cc} = frac.(fields{cc})(:,:,c)';
    end
    
    opts.nrand = 1000;
    opts.minnbchan = 0;
    %opts.parpool = 24; %remove later
    if isempty(plotstats{2,c})
        plotstats{2,c} = EasyClusterCorrect(statdata,datasetinfo,statfun,opts);
    end
    
    
    for cc = 1:length(fields)
        p(2,2,x,y,cc+1).select()
        cluster_topoplot(nanmean(frac.(fields{cc})(:,:,c),1),chanlocs,plotstats{2,c}.prob,plotstats{2,c}.prob < 0.05,1);
        ax(cc) = gca;
        if cc == length(fields)
           pos = ax(cc).Position;
           cbar = colorbar('EastOutside');
           %set(cbar,'Position',[pos(1)+pos(3)-cbar.Position(3) cbar.Position(2:4)]);
        end
    end
    Normalize_Clim(ax,0)
end



AddFigureLabel(p(1,1).axis,'A')
AddFigureLabel(p(2,1,1,1,1).axis,'B')
AddFigureLabel(p(2,2,1,1,1).axis,'C')

