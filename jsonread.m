function [obj] = jsonread(filename)

fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw');
obj = jsondecode(str);