function parsave(filepath,varargin)

for c = 1:2:length(varargin)
    varname = varargin{c};
    input = varargin{c+1};
    eval([varname '=input;']);
end

save(filepath,varname)