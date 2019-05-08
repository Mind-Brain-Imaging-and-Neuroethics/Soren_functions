function [vectvals] = Get_valsforivar(depvar)

varindex = find(contains(outputs.Properties.VariableNames,depvar))

[r,p] = CorrMatrix(allMeasures,outputs,'Plot','off','Threshold','off');

vectvals(1,1) = r(2,varindex);
vectvals(1,2) = p(2,varindex);

[r,p] = CorrMatrix(allMeasuresPCC,outputs,'Plot','off','Threshold','off');

vectvals(2,1) = r(2,varindex);
vectvals(2,2) = p(2,varindex);

[r,p] = CorrMatrix(allMeasuresPACC,outputs,'Plot','off','Threshold','off');

vectvals(3,1) = r(2,varindex);
vectvals(3,2) = p(2,varindex);

[r,p] = CorrMatrix(allMeasuresVC,outputs,'Plot','off','Threshold','off');

vectvals(4,1) = r(2,varindex);
vectvals(4,2) = p(2,varindex);

[r,p] = CorrMatrix(allMeasuresSSC,outputs,'Plot','off','Threshold','off');

vectvals(5,1) = r(2,varindex);
vectvals(5,2) = p(2,varindex);

[r,p] = CorrMatrix(allMeasuresMC,outputs,'Plot','off','Threshold','off');

vectvals(6,1) = r(2,varindex);
vectvals(6,2) = p(2,varindex);