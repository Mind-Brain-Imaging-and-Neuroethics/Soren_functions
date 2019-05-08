function matout = Make_designVect(nbs)

matout = [];
for c = 1:length(nbs)
    matout = [matout ones(1,nbs(c))*c];
end