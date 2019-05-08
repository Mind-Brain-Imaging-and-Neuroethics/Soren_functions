%permutation test thing

rng(100);

seedval = rng;

nperm = 1000;

permP = [];

for q = 1:nperm
    
    shuffleOrder = randperm(31);
    
    allMeasures = allMeasures(shuffleOrder,:);
    selfOutputs = selfOutputs(shuffleOrder,:);
    
    [~,p] = CorrMatrix(allMeasures,selfOutputs,'Plot','off');
    permP = cat(3,permP,p);
end

