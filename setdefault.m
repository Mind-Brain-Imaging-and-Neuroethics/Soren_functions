function [cfg] = setdefault(cfg,field,value)

if ~isfield(cfg,field)
   cfg.(field) = value; 
end