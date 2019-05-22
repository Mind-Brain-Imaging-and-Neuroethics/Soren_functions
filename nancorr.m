function [r,p] = nancorr(in1,in2,varargin)

if CheckInput(varargin,'type')
   type = EasyParse(varargin,'type');
else
    type = 'Pearson';
end

r = zeros(size(in1,2),size(in2,2));
p = r;

for c = 1:size(in1,2)
    for cc = 1:size(in2,2)
        a = in1(:,c);
        a = a(find(~isnan(a)));
        b = in2(:,cc); 
        b = b(find(~isnan(b)));
if ~isempty(a) && ~isempty(b) && (length(a) == length(b)) 
[r(c,cc) p(c,cc)] = corr(vert(a),vert(b),'type',type);
else
r(c,cc) = NaN; 
p(c,cc) = NaN;
end
    end
end

end
