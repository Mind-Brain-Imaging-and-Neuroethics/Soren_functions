function [p] = value_prctile(historicalData,exogenousVariable)

nless = sum(historicalData < exogenousVariable);
nequal = sum(historicalData == exogenousVariable);
p = (nless + 0.5*nequal+1) / length(historicalData);