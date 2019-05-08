%permutation test thing

selfTimeOutputs = selfOutputs(:,5:10);
goodMeasures = allMeasures(:,[1:2 4:5 8:9]);

rng(101);

seedval = rng;

nperm = 5000;

permP = [];

for q = 1:nperm
    
    shuffleOrder = randperm(31);
    
    shuffleMeasures = goodMeasures(shuffleOrder,:);
    %selfOutputs = selfOutputs(shuffleOrder,:);
    
    [~,p] = CorrMatrix(shuffleMeasures,selfTimeOutputs,'Plot','off');
    permP = cat(3,permP,p);
end

for q = 1:nperm
   rankedp(q,:) = sort(reshape(permP(:,:,q),1,[])); 
end

% for q = 1:size(rankedp,2)
%     x = rankedp(:,q);
%     SEM = std(x)/sqrt(length(x));               % Standard Error
%     ts = tinv([0.025  0.975],length(x)-1);      % T-Score
%     CI_rank{q} = mean(x) + ts*SEM;
% end

[~,pvals] = CorrMatrix(goodMeasures,selfTimeOutputs,'Plot','off');
ranked_realp = sort(reshape(pvals,1,[]));

%pvals = reshape(p_corrected,10,10);





