function [PLEest] = MultitaperPLE(onsetSeries,varargin)

range = [exp(-6) exp(0)];

params.Fs = 10;

[S,f] = mtspectrumpt(onsetSeries,params);

f = f';

ind = intersect(find(f > range(1)),find(f < range(2)));

p = polyfit(log(f(ind)),log(S(ind)),1);

PLEest = -p(1);

if length(varargin) > 0 && strcmpi(varargin{1},'Plot')
    y = p(2) + p(1)*log(f(ind));
    plot(log(f(ind)),log(S(ind)));
    hold on;
    plot(log(f(ind)),y,'r--');
    xlabel('Log Frequency')
    ylabel('Log Power')
    title(['Estimated PLE is ' num2str(PLEest)]) 
end