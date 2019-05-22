function [out] = cellstrcmpi(in1,in2)

if iscell(in1)
    for c = 1:length(in1)
        out(c) = strcmpi(in1{c},in2);
    end
else
    for c = 1:length(in2)
        out(c) = strcmpi(in1,in2{c});
    end
end