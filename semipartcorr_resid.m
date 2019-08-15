function [r,p] = semipartcorr_resid(x,y,z)

x = vert(x);
y = vert(y);
if size(z,1) ~= size(y,1)
   z = z'; 
end

[~,~,resid] = regress(y,z);

[r,p] = corr(x,resid);