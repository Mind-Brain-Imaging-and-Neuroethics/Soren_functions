function TTV_ERSP_results_func(settings)


numbands = length(settings.tfparams.fbandnames);

fbands = settings.tfparams.fbandnames;

% Construct allmeas
if ~exist(fullfile(settings.outputdir,[settings.datasetname '_allmeas.mat']),'file') || strcmpi(settings.load_allmeas,'no')
files = dir('*calc.mat');
for c = 1:length(files)
    fprintf([num2str(c) ' '])
    load(files(c).name)
    if c > 1
        dimn = [];
        for q = 1:length(datacalc)
            fields = fieldnames_recurse(datacalc{q});
            fields = cell_unpack(fields);
            
            for cc = 1:length(fields)
                dimn = size(getfield_nest(datacalc{q},fields{cc}));
                %dimn = length(find(size > 1));
                if dimn(end) == 1
                    dimn(end) = [];
                end
                dimn = length(dimn);
                tmp = cat(dimn+1,getfield_nest(allmeas{q},fields{cc}),getfield_nest(datacalc{q},fields{cc}));
                allmeas{q} = assignfield_nest(allmeas{q},fields{cc},tmp);
            end
        end
    else
        allmeas = datacalc;
    end
end
save(fullfile(settings.outputdir,[settings.datasetname '_allmeas.mat']),'allmeas','-v7.3')
elseif strcmpi(settings.load_allmeas,'yes')
load(fullfile(settings.outputdir,[settings.datasetname '_allmeas.mat']))
end

aucindex = settings.aucindex;

alloutputs = struct;

nbchan = settings.nbchan;

%% Comparison of TTV and ERSP time course in each frequency band

for q = 1:numbands
    for c = 1:nbchan
        for i = 1:size(allmeas{q}.ttv.real,3)
            alloutputs.ersp.distrealreal(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.real(c,:,i),[0 1])-NormOntoRange(allmeas{q}.ersp.real(c,:,i),[0 1]));
            alloutputs.ersp.distrealpseudo(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.real(c,:,i),[0 1])-NormOntoRange(allmeas{q}.ersp.pseudo(c,:,i),[0 1]));
            alloutputs.ersp.distpseudoreal(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.pseudo(c,:,i),[0 1])-NormOntoRange(allmeas{q}.ersp.real(c,:,i),[0 1]));
            alloutputs.ersp.distpseudopseudo(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.pseudo(c,:,i),[0 1])-NormOntoRange(allmeas{q}.ersp.pseudo(c,:,i),[0 1]));
            
            alloutputs.itc.distrealreal(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.real(c,:,i),[0 1])-(1-NormOntoRange(allmeas{q}.itc.real(c,:,i),[0 1])));
            alloutputs.itc.distrealpseudo(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.real(c,:,i),[0 1])-(1-NormOntoRange(allmeas{q}.itc.pseudo(c,:,i),[0 1])));
            alloutputs.itc.distpseudoreal(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.pseudo(c,:,i),[0 1])-(1-NormOntoRange(allmeas{q}.itc.real(c,:,i),[0 1])));
            alloutputs.itc.distpseudopseudo(c,i,q) = norm(NormOntoRange(allmeas{q}.ttv.pseudo(c,:,i),[0 1])-(1-NormOntoRange(allmeas{q}.itc.pseudo(c,:,i),[0 1])));
        end
    end
    alloutputs.dist.sigersp(q,:) = signrank_mat(alloutputs.ersp.distrealreal(:,:,q),alloutputs.ersp.distrealpseudo(:,:,q),2);
    alloutputs.dist.sigitc(q,:) = signrank_mat(alloutputs.itc.distrealreal(:,:,q),alloutputs.itc.distrealpseudo(:,:,q),2);
    alloutputs.dist.sigerspvitc(q,:) = signrank_mat(alloutputs.ersp.distrealreal(:,:,q),alloutputs.itc.distrealreal(:,:,q),2);
end


%% Calculation of nonadditivity significance in ERSP

for q = 1:numbands
    alloutputs.ersp.pt.sig(q,:) = signrank_mat(allmeas{q}.naerspindex,zeros(size(allmeas{q}.naerspindex)),2);
    alloutputs.ersp.ttv.sig(q,:) = signrank_mat(allmeas{q}.ttverspindex,zeros(size(allmeas{q}.ttverspindex)),2);
    [r p] = nancorr(allmeas{q}.naerspindex',allmeas{q}.ttverspindex','Type','Spearman');
    alloutputs.ersp.corr.r(:,q) = r(find(eye(nbchan)));
    alloutputs.ersp.corr.p(:,q) = p(find(eye(nbchan)));
    for c = 1:nbchan
try
        alloutputs.ersp.pt.effsize.stats{c,q} = mes(abs((squeeze(trapz(allmeas{q}.naddersp.real(c,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.naddersp.pseudo(c,aucindex,2,:),2))))...
            ,abs((squeeze(trapz(allmeas{q}.naddersp.real(c,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.naddersp.pseudo(c,aucindex,1,:),2)))),'auroc');
        alloutputs.ersp.pt.effsize.vals(c,q) = alloutputs.ersp.pt.effsize.stats{c,q}.auroc;
catch
alloutputs.ersp.pt.effsize.stats{c,q} = NaN;
alloutputs.ersp.pt.effsize.vals(c,q) = NaN;
end
try
        alloutputs.ersp.ttv.effsize.stats{c,q} = mes(allmeas{q}.ttverspindex(c,:),squeeze(trapz(allmeas{q}.ttversp.pseudo(c,aucindex,:),2)),'auroc');
        alloutputs.ersp.ttv.effsize.vals(c,q) = alloutputs.ersp.ttv.effsize.stats{c,q}.auroc;
catch
alloutputs.ersp.ttv.effsize.stats{c,q} = NaN;
alloutputs.ersp.ttv.effsize.vals(c,q) = NaN;
end
    end
end

%% Calculation of ERP nonadditivity

for q = 1:numbands
    alloutputs.erp.pt.sig(q,:) = signrank_mat(allmeas{q}.naerpindex,zeros(size(allmeas{q}.naerpindex)),2);
    alloutputs.erp.ttv.sig(q,:) = signrank_mat(allmeas{q}.ttvindex,zeros(size(allmeas{q}.ttvindex)),2);
    [r p] = nancorr(allmeas{q}.naerpindex',allmeas{q}.ttvindex','Type','Spearman');
    alloutputs.erp.corr.r(:,q) = r(find(eye(nbchan)));
    alloutputs.erp.corr.p(:,q) = p(find(eye(nbchan)));
    for c = 1:nbchan
try
        alloutputs.erp.pt.effsize.stats{c,q} = mes(abs((squeeze(trapz(allmeas{q}.nadderp.real(c,aucindex,2,:),2)) - squeeze(trapz(allmeas{q}.nadderp.pseudo(c,aucindex,2,:),2))))...
            ,abs((squeeze(trapz(allmeas{q}.nadderp.real(c,aucindex,1,:),2)) - squeeze(trapz(allmeas{q}.nadderp.pseudo(c,aucindex,1,:),2)))),'auroc');
        alloutputs.erp.pt.effsize.vals(c,q) = alloutputs.erp.pt.effsize.stats{c,q}.auroc;
catch
alloutputs.erp.pt.effsize.stats{c,q} = NaN;
alloutputs.erp.pt.effsize.vals(c,q) = NaN;
end        
try
alloutputs.erp.ttv.effsize.stats{c,q} = mes(allmeas{q}.ttvindex(c,:),squeeze(trapz(allmeas{q}.ttv.pseudo(c,aucindex,:),2)),'auroc');
        alloutputs.erp.ttv.effsize.vals(c,q) = alloutputs.erp.ttv.effsize.stats{c,q}.auroc;
catch
alloutputs.erp.ttv.effsize.stats{c,q} = NaN;
alloutputs.erp.ttv.effsize.vals(c,q) = NaN;
end   
 end
end

%% Cluster stats
% Do all the cluster stats at the end to minimize time transferring files
% to workers

if ~strcmpi(settings.datatype,'ECoG') || strcmpi(settings.ecog.method,'roi')
    
    ersp_corrstats = cell(1,numbands);
    itc_corrstats = cell(1,numbands);
    
    ersp_diststats = cell(1,numbands);
    itc_diststats = cell(1,numbands);
    
    amp_na_indexcorr_stats = cell(1,numbands);
    amp_na_ttv_stats = cell(1,numbands);
    amp_na_ersp_stats = cell(1,numbands);
    
    ttv_indexstats = cell(1,numbands);
    ttversp_stats = cell(1,numbands);
    ttversp_corrstats = cell(1,numbands);
    
    erp_stats = cell(1,numbands);
    erp_corrstats = cell(1,numbands);
    
    %if isempty(gcp('nocreate'))
    %    parpool(numbands)
    %end
    
    opts = struct;
    opts.nrand = 10000;
    
if strcmpi(settings.datatype,'ECoG')
opts.minnbchan = 0;
end

    for q = 1:numbands
        disp(['Processing ' settings.tfparams.fbandnames{q} ' band'])
        %opts = struct;
        %opts.nrand = 10000;
        
        ttv_indexstats{q} = EasyClusterCorrect({allmeas{q}.ttvindex zeros(size(allmeas{q}.ttvindex))},settings.datasetinfo,'ft_statfun_signrank',opts);
        %ersp_corrstats{q} = EasyClusterCorrect_spearman({allmeas{q}.ttvindex,allmeas{q}.erspindex},settings.datasetinfo);
        %itc_corrstats{q} = EasyClusterCorrect_spearman({allmeas{q}.ttvindex,allmeas{q}.itcindex},settings.datasetinfo);
        
        %ersp_diststats{q} = EasyClusterCorrect_signrank({alloutputs.ersp.distrealreal(:,:,q),alloutputs.ersp.distrealpseudo(:,:,q)},settings.datasetinfo);
        %itc_diststats{q} = EasyClusterCorrect_signrank({alloutputs.itc.distrealreal(:,:,q),alloutputs.itc.distrealpseudo(:,:,q)},settings.datasetinfo);
        ttv_diststats{q} = EasyClusterCorrect({alloutputs.itc.distrealreal(:,:,q),alloutputs.ersp.distrealreal(:,:,q)},settings.datasetinfo,'ft_statfun_signrank',opts);
        
        %amp_na_indexcorr_stats{q} = EasyClusterCorrect_spearman({allmeas{q}.nattvindex.amp,allmeas{q}.naerspindex.amp},settings.datasetinfo);
        %amp_na_ttv_stats{q} = EasyClusterCorrect_signrank({allmeas{q}.nattvindex.amp,zeros(size(allmeas{q}.nattvindex.amp))},settings.datasetinfo);
        amp_na_ersp_stats{q} = EasyClusterCorrect({allmeas{q}.naerspindex,zeros(size(allmeas{q}.naerspindex))},settings.datasetinfo,'ft_statfun_signrank',opts);
        
        %cos_na_indexcorr_stats{q} = EasyClusterCorrect_spearman({allmeas{q}.nattvindex.cos,allmeas{q}.naerspindex.cos},settings.datasetinfo);
        %cos_na_ttv_stats{q} = EasyClusterCorrect_signrank({allmeas{q}.nattvindex.cos,zeros(size(allmeas{q}.nattvindex.cos))},settings.datasetinfo);
        %cos_na_ersp_stats{q} = EasyClusterCorrect_signrank({allmeas{q}.naerspindex.cos,zeros(size(allmeas{q}.naerspindex.cos))},settings.datasetinfo);
        
        ttversp_stats{q} = EasyClusterCorrect({allmeas{q}.ttverspindex zeros(size(allmeas{q}.ttverspindex))},settings.datasetinfo,'ft_statfun_signrank',opts);
        ttversp_corrstats{q} = EasyClusterCorrect({allmeas{q}.naerspindex,allmeas{q}.ttverspindex},settings.datasetinfo,'ft_statfun_spearman',opts);
        
        erp_stats{q} = EasyClusterCorrect({allmeas{q}.naerpindex,zeros(size(allmeas{q}.naerpindex))},settings.datasetinfo,'ft_statfun_signrank',opts);
        erp_corrstats{q} = EasyClusterCorrect({allmeas{q}.naerpindex,allmeas{q}.ttvindex},settings.datasetinfo,'ft_statfun_spearman',opts);
        
    end
    
    %alloutputs.ersp.corrstats = ersp_corrstats;
    %alloutputs.itc.corrstats = itc_corrstats;
    %alloutputs.ersp.diststats = ersp_diststats;
    %alloutputs.itc.diststats = itc_diststats;
    alloutputs.dist.stats = ttv_diststats;
    %alloutputs.amp_na.indexcorr.stats = amp_na_indexcorr_stats;
    %alloutputs.amp_na.ttv.stats = amp_na_ttv_stats;
    alloutputs.ersp.pt.stats = amp_na_ersp_stats;
    % alloutputs.cos_na.indexcorr.stats = cos_na_indexcorr_stats;
    % alloutputs.cos_na.ttv.stats = cos_na_ttv_stats;
    % alloutputs.cos_na.ersp.stats = cos_na_ersp_stats;
    alloutputs.ersp.ttv.stats = ttversp_stats;
    alloutputs.ersp.corr.stats = ttversp_corrstats;
    alloutputs.erp.pt.stats = erp_stats;
    alloutputs.erp.corr.stats = erp_corrstats;
    alloutputs.erp.ttv.stats = ttv_indexstats;
    
    
    
    alloutputs.fdrfields = {'dist.stats','ersp.pt.stats','ersp.ttv.stats','ersp.corr.stats',...
        'erp.pt.stats','erp.ttv.stats','erp.corr.stats'};
elseif strcmpi(settings.datatype,'ECoG') && strcmpi(settings.ecog.method,'mean')
    alloutputs.fdrfields = {'ersp.pt.sig','ersp.ttv.sig','ersp.corr.p','erp.pt.sig','erp.ttv.sig',...
        'erp.corr.p','dist.sigerspvitc'};
end


save([settings.outputdir '/' settings.datasetname '_results.mat'],'alloutputs','-v7.3')

end
