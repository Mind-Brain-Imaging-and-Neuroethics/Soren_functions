function cellout = cell_unpack(cellin)

cellout = [];

for c = 1:length(cellin)
   if ~iscell(cellin{c})
      cellout = [cellout cellin(c)];
   else
       cellout = [cellout cell_unpack(cellin{c})];
   end
end