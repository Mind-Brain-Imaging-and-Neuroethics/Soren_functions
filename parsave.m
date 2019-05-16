function parsave(filepath,varargin)

for c = 1:2:length(varargin)
    varname = varargin{c};
    input = varargin{c+1};
    eval([varname '=input;']);

end
clear input varargin varname c

save(filepath,'-v7.3')
