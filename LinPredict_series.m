function [predictVals] = LinPredict_series(ts,winlen,ple)

rx = CovarFromHurst(ple,ts,var(ts));

for c = (winlen+1):length(ts)
    disp(c)
    [~,predictVals(c)] = LinPredict(ts((c-winlen):(c-1)),rx,ple);
end


disp(c)

