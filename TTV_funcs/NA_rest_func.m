function TTV_ERSP_rest_func(settings)

cd(settings.rest.restdir)
load([settings.outputdir '/' settings.datasetname '_allmeas.mat'])
load([settings.outputdir '/' settings.datasetname '_results.mat'])

%for c = 1:length(settings.tfparams.fbands)
%    if isempty(settings.tfparams.fbands{c})
%        settings.tfparams.fbands{c} = settings.rest.bandpass;
%    end
%end

%eeglab
if isfield(settings.rest,'preprocparams')
    %some preprocessing steps here in EEGLAB
    
    
    cd(settings.rest.preprocdir)
end


%% Calculate bandpower for each frequency band
if strcmpi(settings.datatype,'EEG')
    files = dir(settings.restfiles);
else
    files = dir(settings.restfiles);
end

nbchan = length(settings.datasetinfo.label);
bp = cell(1,length(files));
rel_bp = bp;

restmeas = struct;
filesorder = cell(1,length(files));

rest_fbands = cell(length(files),length(settings.tfparams.fbandnames));
pf =

parfor i = 1:length(files)
    if strcmpi(settings.datatype,'EEG')
        EEG = pop_loadset(files(i).name,pwd);
        data = eeglab2fieldtrip(EEG,'preprocessing','none');
    else
        data = parload(files(i).name,'data');
    end
    
    fbands = settings.tfparams.fbandnames;
    
    if strcmpi(settings.tfparams.pf_adjust,'yes')
        [fbands pf(i)] = NA_convert_alpha_pf(settings,data)
        rest_fbands{i} = horz(fbands);
    end
    
    for c = 1:length(fbands)
        fband = fbands{c};
        if isempty(fband)
            fband = settings.rest.bandpass;
        elseif ~isempty(find(isnan(fband)))
            fband(find(isnan(fband))) = settings.rest.bandpass(1);
        end
        for cc = 1:length(data.label)
            bp{i}(c,cc) = bandpower(data.trial{1}(cc,:),data.fsample,fband);
            rel_bp{i}(c,cc) = bp{i}(c,cc)/bandpower(data.trial{1}(cc,:),data.fsample,settings.rest.bandpass);
        end
    end
    filesorder{i} = files(i).name;
end

if strcmpi(settings.tfparams.pf_adjust,'yes')
    restmeas.rest_fbands = cat(1,rest_fbands{:});
    restmeas.alpha_pf = pf;
end

restmeas.filesorder = filesorder;

bp = cat(3,bp{:});
rel_bp = cat(3,rel_bp{:});

restmeas.bp.vals = bp;
restmeas.rel_bp.vals = rel_bp;

clear bp rel_bp EEG

%% Correlate the resting state bandpower with evoked power,

for q = 1:settings.nfreqs
    [r p] = corr(squeeze(restmeas.bp.vals(q,:,:))',allmeas{q}.erspindex','Type','Spearman');
    restmeas.bp.index.r(:,q) = r(find(eye(nbchan)));
    restmeas.bp.index.p(:,q) = p(find(eye(nbchan)));
    
    [r p] = corr(squeeze(restmeas.bp.vals(q,:,:))',allmeas{q}.naerspindex','Type','Spearman');
    restmeas.bp.naindex.r(:,q) = r(find(eye(nbchan)));
    restmeas.bp.naindex.p(:,q) = p(find(eye(nbchan)));
    
    [r p] = corr(squeeze(restmeas.rel_bp.vals(q,:,:))',allmeas{q}.erspindex','Type','Spearman');
    restmeas.rel_bp.index.r.subject(:,q) = r(find(eye(nbchan)));
    restmeas.rel_bp.index.p.subject(:,q) = p(find(eye(nbchan)));
    
    [r p] = corr(squeeze(restmeas.rel_bp.vals(q,:,:))',allmeas{q}.naerspindex','Type','Spearman');
    restmeas.rel_bp.naindex.r.subject(:,q) = r(find(eye(nbchan)));
    restmeas.rel_bp.naindex.p.subject(:,q) = p(find(eye(nbchan)));
    
    [r p] = corr(mean(squeeze(restmeas.rel_bp.vals(q,:,:)),3),mean(allmeas{q}.erspindex,2),'Type','Spearman');
    restmeas.rel_bp.index.r.electrode(:,q) = r;
    restmeas.rel_bp.index.p.electrode(:,q) = p;
    
    [r p] = corr(mean(squeeze(restmeas.rel_bp.vals(q,:,:)),3),mean(allmeas{q}.naerspindex,2),'Type','Spearman');
    restmeas.rel_bp.naindex.r.electrode(:,q) = r;
    restmeas.rel_bp.naindex.p.electrode(:,q) = p;
    
    restmeas.prestimamp.raw{q} = squeeze(mean(allmeas{q}.raw.ersp(:,settings.real.prestim,:),2));
    restmeas.prestimamp.rel{q} = restmeas.prestimamp.raw{q}./restmeas.prestimamp.raw{1};
    
    [r p] = corr(squeeze(restmeas.bp.vals(q,:,:))',restmeas.prestimamp.raw{q}','Type','Spearman');
    restmeas.bp.prestim.r.subject(:,q) = r(find(eye(nbchan)));
    restmeas.bp.prestim.p.subject(:,q) = p(find(eye(nbchan)));
    
    [r p] = corr(squeeze(restmeas.rel_bp.vals(q,:,:))',restmeas.prestimamp.rel{q}','Type','Spearman');
    restmeas.rel_bp.prestim.r.subject(:,q) = r(find(eye(nbchan)));
    restmeas.rel_bp.prestim.p.subject(:,q) = p(find(eye(nbchan)));
    
    %     [r p] = corr(mean(squeeze(restmeas.bp.vals(q,:,:)),3),mean(restmeas.prestimamp.raw{q},2),'Type','Spearman');
    %     restmeas.bp.prestim.r.electrode(:,q) = r;
    %     restmeas.bp.prestim.p.electrode(:,q) = p;
    %
    %     [r p] = corr(squeeze(restmeas.rel_bp.vals(q,:,:))',mean(restmeas.prestimamp.rel{q},2)','Type','Spearman');
    %     restmeas.rel_bp.prestim.r.electrode(:,q) = r;
    %     restmeas.rel_bp.prestim.p.electrode(:,q) = p;
end


bp_index_stats = cell(1,settings.nfreqs);
bp_naindex_stats = cell(1,settings.nfreqs);
rel_bp_index_stats = cell(1,settings.nfreqs);
rel_bp_naindex_stats = cell(1,settings.nfreqs);
bp_prestim_stats = cell(1,settings.nfreqs);
rel_bp_prestim_stats = cell(1,settings.nfreqs);
%bp_mediation_raw = cell(1,settings.nfreqs);
%bp_mediation_cluster = cell(1,settings.nfreqs);
bp_mediation_stats = cell(1,settings.nfreqs);
rel_bp_mediation_stats = cell(1,settings.nfreqs);


opts = struct;
opts.display = 0;
opts.verbose = 0;
opts2 = struct;
opts2.nrand = 10000;

parfor q = 1:settings.nfreqs
    rel_bp_index_stats{q} = EasyClusterCorrect_spearman({squeeze(restmeas.rel_bp.vals(q,:,:)),allmeas{q}.erspindex},settings.datasetinfo,opts2);
    rel_bp_naindex_stats{q} = EasyClusterCorrect_spearman({squeeze(restmeas.rel_bp.vals(q,:,:)),allmeas{q}.naerspindex},settings.datasetinfo,opts2);
    
    rel_bp_prestim_stats{q} = EasyClusterCorrect_spearman({squeeze(restmeas.rel_bp.vals(q,:,:)),restmeas.prestimamp.rel{q}},settings.datasetinfo,opts2);
    
    if isfield(rel_bp_prestim_stats{q},'posclusters') && ~isempty(rel_bp_prestim_stats{q}.posclusters) && ~isempty(find(extractfield(rel_bp_prestim_stats{q}.posclusters,'prob') < 0.05))
        rel_bp_mediation_stats{q} = mediationAnalysis0(double(mean(allmeas{q}.erspindex(find(rel_bp_prestim_stats{q}.mask),:),1))',...
            double(squeeze(mean(restmeas.rel_bp.vals(q,find(rel_bp_prestim_stats{q}.mask),:),2))),...
            double(mean(restmeas.prestimamp.rel{q}(find(rel_bp_prestim_stats{q}.mask),:),1))',opts);
    end
    if (isfield(rel_bp_prestim_stats{q},'negclusters') && ~isempty(rel_bp_prestim_stats{q}.negclusters) && ~isempty(find(extractfield(rel_bp_prestim_stats{q}.negclusters,'prob') < 0.05)))
        rel_bp_mediation_stats{q} = mediationAnalysis0(double(mean(allmeas{q}.erspindex(find(rel_bp_prestim_stats{q}.mask),:),1))',...
            double(squeeze(mean(restmeas.rel_bp.vals(q,find(rel_bp_prestim_stats{q}.mask),:),2))),...
            double(mean(restmeas.prestimamp.rel{q}(find(rel_bp_prestim_stats{q}.mask),:),1))',opts);
    end
end


restmeas.rel_bp.index.stats = rel_bp_index_stats;
restmeas.rel_bp.naindex.stats = rel_bp_naindex_stats;
restmeas.rel_bp.prestim.stats = rel_bp_prestim_stats;
restmeas.rel_bp.mediation = rel_bp_mediation_stats;

restmeas.fdrfields = {'rel_bp.index.stats','rel_bp.naindex.stats','rel_bp.prestim.stats','rel_bp.mediation'};

save([settings.outputdir '/' settings.datasetname '_restmeas.mat'],'restmeas')


end
