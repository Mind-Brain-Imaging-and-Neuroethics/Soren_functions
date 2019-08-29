function pvals = Permtest_anovan(y,group,nrand)


[~,tmain] = anovan(y,group,'model','interaction','display','off');

ngroups = length(group);

fmain = cat(1,tmain{2:(1+sum(1:ngroups)),6});

for i = 1:nrand
    perm = randperm(length(y));
    for c = 1:length(group)
        permgroup{c} = group{c}(perm);
    end
    [~,tperm] = anovan(y,permgroup,'model','interaction','display','off');
    fperm(:,i) = cat(1,tperm{2:(1+sum(1:ngroups)),6});
end

for c = 1:size(fmain,1)
    pvals(c) = value_prctile(fperm(c,:),fmain(c));
    if pvals(c) > 0.5
       pvals(c) = 1-pvals(c); 
    end
    pvals(c) = pvals(c)*2; %two tailed test
end