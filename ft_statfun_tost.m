function s = ft_statfun_tost(cfg,dat,design)

tmp = unique(design);
indices{1} = find(design == tmp(1));
indices{2} = find(design == tmp(2));

p = ones(size(dat,1),1);
for c = 1:size(dat,1)
    [p1,p2] = TOST(dat(c,indices{1}),dat(c,indices{2}),cfg.critval(1),cfg.critval(2));
    teststat(c) = max(p1,p2);
    if mean(dat(c,indices{1})) > mean(dat(c,indices{2}))
        teststat(c) = 1-teststat(c);
    else
        teststat(c) = -1+teststat(c);
    end
end

s.critval = [-0.95 0.95];
