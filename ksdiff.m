function [D] = ksdiff(x1,x2)

% Calculates the Kolmogorov-Smirnov test statistic for two vectors
% representing probability density functions

cx1 = cumsum(x1)/sum(x1);
cx2 = cumsum(x2)/sum(x2);

cx1 = cx1(1:end-1);
cx2 = cx2(1:end-1);

deltaCDF = abs(cx1 - cx2);

D = max(deltaCDF);


