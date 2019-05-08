function [p] = Diff_Var(x,y)

xdiff = Get_AllDifs(x);
ydiff = Get_AllDifs(y);

p = ranksum(xdiff,ydiff);