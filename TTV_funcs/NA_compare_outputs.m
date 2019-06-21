function compstats = NA_compare_outputs(settings,allmeas1,allmeas2)

opts.nrand = 1000; 
opts.minnbchan = 1;
opts.parpool = settings.pool;

for c = 1:settings.nfreqs
    pt_diff_stats{c} = EasyClusterCorrect({permute(squeeze(allmeas1{c}.naddersp.diff(:,:,2,:)-allmeas1{c}.naddersp.diff(:,:,1,:)),[1 3 2]),...
                permute(squeeze(allmeas2{c}.naddersp.diff(:,:,2,:)-allmeas2{c}.naddersp.diff(:,:,1,:)),[1 3 2])},...
            settings.datasetinfo,'ft_statfun_fast_signrank',opts);
    ttv_diff_stats{c} = EasyClusterCorrect({permute(allmeas1{c}.ttversp.real,[1 3 2]),permute(allmeas2{c}.ttversp.real,[1 3 2])},...
         settings.datasetinfo,'ft_statfun_fast_signrank',opts);
     
    ersp_diff_stats{c} = EasyClusterCorrect({permute(allmeas1{c}.ersp.real,[1 3 2]) permute(allmeas2{c}.ersp.real,[1 3 2])},...
        settings.datasetinfo,'ft_statfun_fast_signrank',opts);
end

compstats.ptdiff = pt_diff_stats;
compstats.ttvdiff = ttv_diff_stats;
compstats.erspdiff = ersp_diff_stats;