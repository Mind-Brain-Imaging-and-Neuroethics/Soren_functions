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
cfg.toi = [-1:0.02:1.5]; cfg.parflag = 'yes';

parpool(22)

for i = indices
    if ~exist([files(i).name '_IRASAtf.mat'],'file')
    load(files(i).name)
    [osci,specs] = IRASA_tf(cfg,data);
    cfg.oscifrac = 'frac';
    frac = IRASA_tf(cfg,data,specs);
    save([files(i).name '_IRASAtf.mat'],'osci','frac','specs')
    end
end

end
