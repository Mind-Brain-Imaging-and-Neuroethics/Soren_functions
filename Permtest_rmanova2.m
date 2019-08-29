function pvals = Permtest_rmanova2(y,subject,group1,group2,names,nrand)

tmain = rm_anova2(y,subject,group1,group2,names);

ngroups = 2;

fmain = cat(1,tmain{2:(1+sum(1:ngroups)),5});

for i = 1:nrand
    perm = randperm(length(group1));
    permgroup1 = group1(perm);
    permgroup2 = group2(perm);
    permsub = subject(perm);
    tperm = rm_anova2(y,permsub,permgroup1,permgroup2,names);
    fperm(:,i) = cat(1,tperm{2:(1+sum(1:ngroups)),5});
end

for c = 1:size(fmain,1)
    pvals(c) = value_prctile(fperm(c,:),fmain(c));
    if pvals(c) > 0.5
       pvals(c) = 1-pvals(c); 
    end
    pvals(c) = pvals(c)*2; %two tailed test
end