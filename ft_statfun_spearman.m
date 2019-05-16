function [s, cfg] = ft_statfun_spearman(cfg, dat, design)

tmp = unique(design);
indices{1} = find(design == tmp(1));
indices{2} = find(design == tmp(2));

p = ones(size(dat,1),1);
for c = 1:size(dat,1)
    tic
    [r,p(c)] = nancorr(dat(c,indices{1})',dat(c,indices{2})','Type','Spearman');
    if r > 0
        p(c) = 1-p(c);
    else
        p(c) = -1+p(c);
    end
    toc
end

s.stat = p;

switch cfg.tail
    case 0
        s.critval = [-0.95 0.95];
    case 1
        s.critval = 0.95;
    case -1
        s.critval = -0.95;
end
s.df = size(dat,2)-2;
end
