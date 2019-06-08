function TTV_ERSP_behaviour_func(settings)

load([settings.outputdir '/' settings.datasetname '_allmeas.mat'])
load([settings.outputdir '/' settings.datasetname '_results.mat'])
if isfield(settings,'rest')
    load([settings.outputdir '/' settings.datasetname '_restmeas.mat'])
end

prestim_pseudo = settings.pseudo.prestim;
prestim_real = settings.real.prestim;
poststim_pseudo = settings.pseudo.poststim;
poststim_real = settings.real.poststim;

for i = 1:length(settings.behav)
    %% Load behavioural data and set things up
    disp(['Performing behavioural analysis ' num2str(i) '...'])
    if settings.behav{i}.file(1) == '/'
        load(settings.behav{i}.file)
    else
        load([settings.outputdir '/' settings.behav{i}.file])
    end
    
    if iscell(behav_data)
        settings.behav{i}.design = 'within';
    else
        settings.behav{i}.design  = 'between';
    end
    indvar = [];
    foi = Subset_index(settings.tfparams.fbandnames,settings.behav{i}.foi);
    %% Get the independent variable
    
    if exist([settings.outputdir '/indvar' num2str(i) '.mat'],'file') && isfield(settings.behav{i},'loadvar') && strcmpi(settings.behav{i}.loadvar,'yes')
        indvar = parload([settings.outputdir '/indvar' num2str(i) '.mat'],'indvar');
    else
        switch settings.behav{i}.indvar
            case 'erspindex'
                if strcmpi(settings.behav{i}.design,'within')
                    files = dir(fullfile(settings.inputdir,settings.files));
                    
                    indvar = cell(1,length(files));
                    parfor q = 1:length(files)
                        if strcmpi(settings.datatype,'EEG')
                            EEG = pop_loadset(files(q).name,files(q).folder);
                            data = eeglab2fieldtrip(EEG,'preprocessing','none');
                            EEG = [];
                        else
                            data = parload(fullfile(files(q).folder,files(q).name),'data');
                        end
                        
                        timefreq_data = [];
                        tmpsettings = settings;
                        if strcmpi(settings.tfparams.pf_adjust,'yes')
                            freqs = settings.tfparams.fbands(q,foi);
                        else
                            freqs = settings.tfparams.fbands(foi);
                        end
                        %tmpsettings.tfparams.fbands = settings.tfparams.fbands(foi);
                        tmpsettings.tfparams.fbandnames = settings.tfparams.fbandnames(foi);
                        timefreq_data = NA_get_tfdata(tmpsettings,data,freqs);
                        for qq = 1:length(foi)
                            for c = 1:length(timefreq_data{qq}.trial)
                                ersp_real = abs(timefreq_data{qq}.trial{c}(:,poststim_real)) - ...
                                    mean(abs(timefreq_data{qq}.trial{c}(:,prestim_real)),2);
                                %ersp_pseudo = abs(timefreq_data{foi(qq)}.trial{c}(:,poststim_pseudo)) - ...
                                %   mean(abs(timefreq_data{foi(qq)}.trial{c}(:,prestim_pseudo)),2);
                                switch settings.units
                                    case 'prcchange'
                                        %  ersp_pseudo = 100*ersp_pseudo./mean(abs(timefreq_data{foi(qq)}.trial{c}(:,prestim_pseudo)),2);
                                        ersp_real = 100*ersp_real./mean(abs(timefreq_data{qq}.trial{c}(:,prestim_real)),2);
                                    case 'zscore'
                                        % ersp_pseudo = zscore(ersp_pseudo,0,2);
                                        ersp_real = zscore(ersp_real,0,2);
                                    case 'log'
                                        %ersp_pseudo = 10*log10(ersp_pseudo);
                                        ersp_real = 10*log10(ersp_real);
                                end
                                indvar{q}(c,:,qq) = trapz(ersp_real,2);
                                %indvar{q}(c,:,qq) = trapz(ersp_real-ersp_pseudo,2);
                            end
                        end
                    end
                else
                    for q = 1:size(allmeas{1}.erspindex,2)
                        for qq = 1:length(foi)
                            indvar(:,q,qq) = allmeas{foi(qq)}.erspindex(:,q);
                        end
                    end
                end
                
            case 'ttvindex'
                for q = 1:size(allmeas{1}.ttvindex,2)
                    for qq = 1:length(foi)
                        indvar(:,q,qq) = allmeas{foi(qq)}.ttvindex(:,q);
                    end
                end
                
                
            case 'rest_bp'
                for q = 1:size(allmeas{1}.erspindex,2)
                    for qq = 1:length(foi)
                        indvar(:,q,qq) = restmeas.bp.vals(foi(qq),:,q);
                    end
                end
            case 'rest_relbp'
                for q = 1:size(allmeas{1}.erspindex,2)
                    for qq = 1:length(foi)
                        indvar(:,q,qq) = restmeas.rel_bp.vals(foi(qq),:,q);
                    end
                end
        end
    end
    save([settings.outputdir '/indvar' num2str(i) '.mat'],'indvar')
    %% Do the analysis
    if strcmpi(settings.behav{i}.design,'within')
        if islogical(behav_data{1}) || (length(unique(behav_data{1})) == 2)
            settings.behav{i}.datatype = 'categorical';
        else
            settings.behav{i}.datatype = 'continuous';
        end
        % delete(gcp('nocreate'))
        % parpool
        for q = 1:length(foi)
            allindvar = [];
            nbtrials = cell2mat(cellfun(@(datain)size(datain,1),indvar,'UniformOutput',false));
            subs = Make_designVect(nbtrials);
            for qq = 1:settings.nbchan
                indvar_vect = [];
                for c = 1:length(indvar)
                    indvar_vect = [indvar_vect indvar{c}(:,qq,q)'];
                end
                allindvar(qq,:) = indvar_vect;
                
                depvar_vect = [];
                for c = 1:length(behav_data)
                    if size(behav_data{c},1) > 1
                        depvar_vect = [depvar_vect behav_data{c}'];
                    else
                        depvar_vect = [depvar_vect behav_data{c}];
                    end
                end
                
                if isfield(settings.behav{i},'threshold')
                    badindices = [find(depvar_vect < settings.behav{i}.threshold(1)) find(depvar_vect > settings.behav{i}.threshold(2))];
                end
                
                designMat = [subs' indvar_vect' depvar_vect'];
                designTbl = array2table(double(designMat),'VariableNames',{'Subject',settings.behav{i}.indvar,'Behav'});
                if strcmpi(settings.behav{i}.datatype,'categorical')
                    model{qq} = fitglme(designTbl,['Behav ~ ' settings.behav{i}.indvar ...
                        ' + (' settings.behav{i}.indvar '|Subject)'],'Distribution','binomial','Link','logit');
                else
                    model{qq} = fitlme(designTbl,['Behav ~ ' settings.behav{i}.indvar ...
                        ' + (' settings.behav{i}.indvar '|Subject)']);
                end
            end
            behav_meas{i}.model{q} = model;
            
            depvar_vect = [];
            for c = 1:length(behav_data)
                if size(behav_data{c},1) > 1
                    depvar_vect = [depvar_vect behav_data{c}'];
                else
                    depvar_vect = [depvar_vect behav_data{c}];
                end
            end
            %behav_meas{i}.cluster{q} = EasyClusterCorrect_mixedeff({allindvar repmat(depvar_vect,settings.nbchan,1)},...
            %    subs,settings.datasetinfo,'dep~ind+(ind|subs)');
        end
        
    elseif strcmpi(settings.behav{i}.design,'between')
        for q = 1:length(foi)
            if size(behav_data,1) == 1
                behav_data = behav_data';
            end
            [r p] = corr(indvar(:,:,q)',behav_data,'Type','Spearman');
            behav_meas{i}.r(:,q) = r;
            behav_meas{i}.p(:,q) = p;
            opts.nrand = 10000;
            opts.minnbchan = 1;
            behav_meas{i}.stats = EasyClusterCorrect({indvar(:,:,q) repmat(behav_data,1,settings.nbchan)'},settings.datasetinfo,'ft_statfun_correlationT',opts);
        end
    end
    behav_meas{i}.settings = settings.behav{i};
    behav_meas{i}.indvar = indvar;
    behav_meas{i}.depvar = behav_data;
end


save([settings.outputdir '/' settings.datasetname '_behav.mat'],'behav_meas')

end
