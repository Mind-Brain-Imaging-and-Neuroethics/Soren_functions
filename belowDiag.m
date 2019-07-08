function outmat = belowDiag(inmat)

mask = tril(ones(size(inmat)))-eye(size(inmat));
outmat = inmat.*mask;