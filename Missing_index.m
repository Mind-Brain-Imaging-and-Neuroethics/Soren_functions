function [missing_indx] = Missing_index(cell_all,cell_sub)

% count = 1; 
% for c = 1:length(cell_all)
%    if ~isempty(find(strcmpi(cell_all{c},cell_sub)))
%        indx(count) = c; 
%        count = count+1;
%    end
% end

sub_indx = Subset_index(cell_all,cell_sub);
all_indx = 1:length(cell_all);
missing_indx = find(~ismember(all_indx,sub_indx));
