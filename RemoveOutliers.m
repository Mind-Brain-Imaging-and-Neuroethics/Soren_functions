function [c,d,rmout] = RemoveOutliers(a,b)

rmout = [find(isoutlier(a))' find(isoutlier(b))'];
a(rmout) = [];
b(rmout) = [];
c = a;
d = b;