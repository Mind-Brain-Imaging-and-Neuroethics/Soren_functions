function [output] = getdim(A,dim)

sz = size(A);
inds = repmat({1},1,ndims(A));
inds{dim} = 1:sz(dim);
output = A(inds{:});