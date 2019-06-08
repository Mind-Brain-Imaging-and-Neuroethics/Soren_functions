function settings = NA_alpha_pf(settings,filesuffix)

files = dir(files);

pool = gcp('nocreate');

if isempty(pool) || pool.NumWorkers ~= settings.pool
    parpool(settings.pool)
end

fbands_alpha = cell(1,length(files));
fbands_theta = cell(1,length(files));
fbands_beta = cell(1,length(files));
fbands = repmat(horz(settings.tfparams.fbands),length(files),1);

individ = {'Theta','Alpha','Beta'};
individindx = Subset_index(settings.tfparams.fbandnames,individ);

parfor i = 1:length(files)
    
    if strcmpi(settings.datatype,'EEG')
        EEG = pop_loadset('filename',files(i).name,'filepath',files(i).folder);
        data = eeglab2fieldtrip(EEG,'preprocessing','none');
    else
        data = parload(files(i).name,'data');
        data = ft_concat(data);
    end
    
    [psum,~,f] = restingIAF(data.trial{1},length(data.label),3,[1 40],500,[7 14],11,5);
    
    if ~isnan(psum.paf)
        pfs(i) = psum.paf
    else
        pfs(i) = psum.cog;
    end
    
    fbands_alpha{i} = [f(psum.iaw(1)) f(psum.iaw(2))];
    fbands_theta{i} = [4 f(psum.iaw(1))];
    fbands_beta{i} = [f(psum.iaw(2)) 30];
end

fbands(:,individindx(1)) = vert(fbands_theta);
fbands(:,individindx(2)) = vert(fbands_alpha);
fbands(:,individindx(3)) = vert(fbands_beta);

settings.tfparams.fbands = fbands;
settings.tfparams.alpha_pf = pfs;
