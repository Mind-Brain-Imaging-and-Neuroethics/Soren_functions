function [PV] = Prop_Var(datain)

alldifs = Get_AllDifs(datain);

PV = sum(alldifs)/length(alldifs);