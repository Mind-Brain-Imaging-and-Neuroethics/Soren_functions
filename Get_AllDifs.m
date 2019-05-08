function [alldifs,alldifsmat] = Get_AllDifs(vectin)

for c = 1:length(vectin)
    for cc = 1:c
        if max(vectin(c),vectin(cc))
            alldifsmat(c,cc) = abs(vectin(c) - vectin(cc))/abs(max(vectin(c),vectin(cc)));
        else
            alldifsmat(c,cc) = 0;
        end
    end
end

mask = tril(ones(size(alldifsmat))) - eye(size(alldifsmat));
alldifs = reshape(alldifsmat(find(mask)),[],1);

end