function [output] = CovarFromHurst(PLE, timeSeries, V,varargin)

K = length(timeSeries);
if length(varargin) < 1
    output = zeros(1,K);
else
    value = cell2mat(varargin(1));
end


%V = var(timeSeries);
hExp = (PLE+1)/2;
if length(varargin) > 0
    output = (V/2)*((abs(value+1)^(2*hExp))-(2*abs(value)^(2*hExp))...
        +(abs(value-1)^(2*hExp)));
else
    for k = 0:K
        output(k+1) = (V/2)*(abs(k+1)^(2*hExp)-2*abs(k)^(2*hExp)...
            +abs(k-1)^(2*hExp));
    end
end

