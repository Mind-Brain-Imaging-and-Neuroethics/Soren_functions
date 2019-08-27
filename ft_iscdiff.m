function [stats,tbl] = ft_iscdiff(cfg,data)

if isempty(cfg)
    cfg = struct;
end

cfg = setdefault(cfg,'cond',{'Condition1' 'Condition2'});
cfg = setdefault(cfg,'thresh',0.1);

for c = 1:length(data{1}.meas)
   
   p(c) = Permtest_ISC(data{1}.data(:,:,c)',data{2}.data(:,:,c)',1000,'spearman');
   isc1{c} = corr(data{1}.data(:,:,c)','Type','Spearman');
   isc2{c} = corr(data{2}.data(:,:,c)','Type','Spearman');
   meanisc1(c) = mean(mean(isc1{c}));
   meanisc2(c) = mean(mean(isc2{c}));
   diff(c) = meanisc1(c)-meanisc2(c);
   if diff(c) > 0
      dirn{c} = [cfg.cond{1} ' > ' cfg.cond{2}]; 
   else
       dirn{c} = [cfg.cond{2} ' > ' cfg.cond{1}];
   end
end

stats = struct;
stats.p = p;
stats.isc1 = isc1;
stats.isc2 = isc2;
stats.meas = data{1}.meas;

tbl = table;
pmask = find(p < cfg.thresh);
tbl.meas = vert(cellfun(@func2str,data{1}.meas(pmask),'UniformOutput',false));
tbl.pvalue = vert(p(pmask));
tbl.([cfg.cond{1} '_mean']) = vert(meanisc1(pmask));
tbl.([cfg.cond{2} '_mean']) = vert(meanisc2(pmask));
tbl.diff = vert(diff(pmask));
tbl.summary = vert(dirn(pmask));