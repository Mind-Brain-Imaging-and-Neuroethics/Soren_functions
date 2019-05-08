function [outliers] = isoutlier2(data,threshold)

data = zscore(data);

outliers = data > threshold;