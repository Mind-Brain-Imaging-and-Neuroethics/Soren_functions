function s = ft_statfun_tost(cfg,dat,design)

tmp = unique(design);
indices{1} = find(design == tmp(1));
indices{2} = find(design == tmp(2));

teststat = ones(size(dat,1),1);
for c = 1:size(dat,1)
    [p1,p2] = TOST(dat(c,indices{1}),dat(c,indices{2}),cfg.eqinterval(1),cfg.eqinterval(2),cfg.alpha);
    teststat(c) = max(p1,p2);
    if mean(dat(c,indices{1})) > mean(dat(c,indices{2}))
        teststat(c) = 1-teststat(c);
    else
        teststat(c) = -1+teststat(c);
    end
end

s.stat = teststat;

s.critval = [-0.95 0.95];

s.df = size(dat,2)-2;
