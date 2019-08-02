function cellout = catcell(varargin)

cellout = varargin{1};

for c = 2:length(varargin)
   for cc = 1:length(varargin{c})
      cellout{cc} = [cellout{cc} varargin{c}{cc}]; 
   end
end