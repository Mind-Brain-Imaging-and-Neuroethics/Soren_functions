function CamCAN_IRASA_tf(indx)
% script for getting IRASA time-frequency stuff

addpath(genpath('/project/def-gnorthof/sorenwt/MATLAB'))
addpath('/project/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

cd /scratch/sorenwt/camcan/Preprocessed/Task/Epoched/

files = dir('*1Hz.mat');

indices = ((indx-1)*32+1):(indx*32);

indices(indices > length(files)) = [];

cfg.oscifrac = 'osci'; cfg.winsize = 1.5; cfg.modifywindow = 'no';
cfg.toi = [-1.25:0.02:0.75]; cfg.parflag = 'yes';

pc = parcluster('local');

pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_ARRAY_TASK_ID'));

parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

for i = indices
    if ~exist([files(i).name '_IRASAtf.mat'],'file')
        load(files(i).name)
        
        [osci,specs] = IRASA_tf(cfg,data);
        save([files(i).name '_IRASAtf.mat'],'osci','specs','-v7.3');
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

end
