function settings = NA_get_alpha_pf(settings)

pool = gcp('nocreate');

if isempty(pool) || pool.NumWorkers ~= settings.pool
    delete(gcp('nocreate'))
    parpool(settings.pool)
end

files = dir(fullfile(settings.inputdir,settings.files));

fbands = cell(length(files),1);

pf = zeros(length(files),1);

parfor i = 1:length(files)
      if strcmpi(settings.datatype,'EEG')
        EEG = pop_loadset('filename',files(i).name,'filepath',files(i).folder);
        data = eeglab2fieldtrip(EEG,'preprocessing','none');
    else
        data = parload(files(i).name,'data');
        data = ft_concat(data);
      end
    
      [band tmppf] = NA_convert_alpha_pf(settings,data);
      fbands{i} = horz(band);
      pf(i) = tmppf;

end

fbands = cat(1,fbands{:});
settings.tfparams.fbands = fbands;
settings.alpha_pf = pf;


