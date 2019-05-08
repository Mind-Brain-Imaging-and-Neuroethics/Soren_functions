function ft_topoplot_vec(layout,vect,label)

tlock = [];
tlock.avg = vect;
tlock.dimord = 'chan_time';
tlock.label = label;
tlock.time = 1;

if ischar(layout)
    cfg = []; cfg.layout = layout;
    layout = ft_prepare_layout(cfg);
end
cfg = [];
cfg.layout = layout;
cfg.interpolateovernan = 'yes';
cfg.lay = layout;
cfg.channel = label;

ft_topoplotER(cfg,tlock);