addpath(genpath('/home/soren/Documents/MATLAB'))
addpath('/group/northoff/share/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

cd /home/soren/Documents/camcan/Preprocessed/Task/Epoched/

files = dir('*epoched.mat');

indices = [29:32:length(files) 30:32:length(files)];

indices(indices > length(files)) = [];

cfg.oscifrac = 'osci'; cfg.winsize = 1.5; cfg.modifywindow = 'no';
cfg.toi = [-1.25:0.02:0.75]; cfg.parflag = 'yes';

%pc = parcluster('local');

%pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_ARRAY_TASK_ID'));

%parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

for i = indices
    disp(['Processing ' files(i).name '...'])
    if ~exist([files(i).name '_IRASAtf.mat'],'file')
        load(files(i).name)
delete(gcp('nocreate'))
parpool(24)


%	for c = 1:length(data.time)
%	     data.time{c} = linspace(-2,1.5,length(data.time{c}));
%	end
        
        [osci,specs] = IRASA_tf(cfg,data);
        save([files(i).name '_IRASAtf.mat'],'osci','-v7.3');
        osci = [];

        cfg.oscifrac = 'frac';
        frac = IRASA_tf(cfg,data,specs);
        save([files(i).name '_IRASAtf.mat'],'frac','-v7.3','-append');
        frac = [];
        
        cfg.oscifrac = 'mixd';
        mixd = IRASA_tf(cfg,data,specs);
        save([files(i).name '_IRASAtf.mat'],'mixd','-v7.3','-append')
        
        specs = [];
        mixd = [];
    end
end
