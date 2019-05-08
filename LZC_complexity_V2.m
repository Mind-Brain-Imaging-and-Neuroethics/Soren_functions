function [LZComplexity] = LZC_compression(series)
% LZC calculates the Lempel-Ziv complexity of a sequence.
%
% Version: 2.1
%
% Creation:   18th April 2006
% Last modification: 22th February 2015
%


%make sure the series is a column vector

if size(series,2) > 1
    series = series';
end

%convert series to binary representation

minV = min(series);
maxV = max(series);

edges = linspace(minV,maxV,16);

series = discretize(series,edges);

%series = de2bi(series);

a = [];
for i = 1:length(series)
   a = horzcat(a,series(i,:));
end

series = a;

%now set the parameters
n = length(series);
b = n/log2(n);

%% LZ algorithm
% Initialization
c = 1;	% Complexity is set to 1
S = series(1);
Q = series(2);

for i = 2:n
    SQ = [S,Q];
	SQ_pi = [SQ(1:(length(SQ)-1))];
    
    % Look for pattern Q inside SQ_pi
	k = findstr(Q,SQ_pi); 

	if length(k) == 0
   	    % New sequence
        % New value of complexity
   	    c = c+1;				
        if (i+1)>n
            break;
        else
            % Updating
            S = [S,Q];
            Q = series(i+1);
        end
   else
        % Not new sequence
        if (i+1)>n
            break;
        else
            % Updating
            Q = [Q,series(i+1)];	
   	    end
    end
end

% Normalization (LZComplexity must be independent of the length of the series)
LZComplexity = c/b;
