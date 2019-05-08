function [PLEest] = SimplePLE(timeSeries,varargin)

range = [0.1 70];

[Pxx,F] = pwelch(timeSeries);

ind = intersect(find(F > range(1)),find(F < range(2)));

p = polyfit(log(F(ind)),log(Pxx(ind)),1);

PLEest = -p(1);

if length(varargin) > 0 && any(varargin{1} == 'Plot')
    y = p(2) + p(1)*log(F(ind));
    plot(log(F(ind)),log(Pxx(ind)));
    hold on;
    plot(log(F(ind)),y,'r--');
    xlabel('Log Frequency')
    ylabel('Log Power')
    title(['Estimated PLE is ' num2str(PLEest)]) 
end
