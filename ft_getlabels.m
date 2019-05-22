function labels = ft_getlabels(data,atlas)

for c = 1:length(data.label)
    cfg = []; cfg.roi = data.elec.elecpos(c,:); cfg.atlas = atlas; cfg.inputcoord = atlas.coordsys; cfg.output = 'single';
    cfg.maxqueryrange = 3;
    try
        tmp = ft_volumelookup(cfg,atlas);
        labels{c} = tmp.name(find(tmp.count == 1));
        if isempty(labels{c})
            labels{c} = 'no_label_found';
        end
    catch
        labels{c} = 'no_label_found';
    end
end
