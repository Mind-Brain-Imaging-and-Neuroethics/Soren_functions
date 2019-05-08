function [indx] = Subset_index(cell_all,cell_sub)

indx = zeros(1,length(cell_sub));

for c = 1:length(cell_sub)
   indx(c) = find(strcmpi(cell_sub{c},cell_all)); 
end