function [CI_lo,CI_hi]=ci_powertail(data,alpha,tail,clt,plotflag)
%CI_POWERTAIL Test for 'dragon kings' vs. 'black swans'.
%   [CI_LO,CI_HI]=CI_POWERTAIL(DATA,ALPHA,TAIL) calculates confidence 
%   intervals of significance level ALPHA for a power law fitted to the 
%   right tail of the empirical distrubution function (edf) of the 
%   observations in DATA. Only observations lying in the interval 
%   [TAIL(1)-quantile of DATA,TAIL(2)-quantile of DATA] are used to obtain 
%   the least squares fit [1].
%   CI_POWERTAIL(DATA,ALPHA,TAIL,CLT) allows to choose between the exact
%   but slower to compute (binomial law based; CLT=0) and asymptotic 
%   (Gaussian law based; CLT=1, default) confidence intervals. 
%   CI_POWERTAIL(DATA,ALPHA,TAIL,CLT,1) plots the power law and the 
%   intervals in the current figure. Observations lying outside the curves 
%   spanned by the confidence intervals are likely to be (i.e. with 
%   probability 1-ALPHA) 'dragon kings', i.e. deviations from power law 
%   events (commonly known as 'black swans'). Note, that Taleb [2] calls 
%   these 'tractable scientifically' events 'Mandelbrotian gray swans', to 
%   distinguish them from the 'totally intractable' black swans.
%
%   Sample use:
%     >> x = 2*(rand(5000,1).^(-1/1) - 1); 
%     >> [ci_lo,ci_hi] = ci_powertail(x,.05,[.1 .01],1,1);
%
%   References:
%   [1] J. Janczura, R. Weron (2012) Black swans or dragon kings? A simple 
%       test for deviations from the power law, European Physical Journal 
%       - Special Topics (EPJ ST) 205, 79-93 
%       (http://dx.doi.org/10.1140/epjst/e2012-01563-9). Working paper  
%       available from: http://ideas.repec.org/p/wuu/wpaper/hsc1101.html
%   [2] N.N.Taleb (2007) The Black Swan: The Impact of the Highly 
%       Improbable, Random House.

%   Written by Joanna Janczura and Rafal Weron (2012.03.03). 
%   Based on functions powerfit.m, ci_powerlaw.m, ci_powerlaw_clt.m
%   written by Joanna Janczura, Rafal Weron (1999.09.09, 2011.04.26). 

if nargin<5
    plotflag = 0;
end
if nargin<4
    clt = 1;
end
n = length(data);
k = ceil(tail(1)*n);
s = ceil(tail(2)*n);

% Compute edf
[y,x] = ecdf(data);
% Shift edf by 1/2 of the step
y = y - 1/(2*n);

% Fit a power law to the tail of the edf
[power,coeff] = powerfit(x(end-k:end-s),1 - y(end-k:end-s));
powertail = coeff*x(end-k:end).^power;
if clt, % Asymptotic, Gaussian CI
    CI_lo = max(0, powertail - sqrt(powertail.*(1 - powertail)./n).*norminv(1 - alpha/2));
    CI_hi = powertail + sqrt(powertail.*(1 - powertail)./n).*norminv(1 - alpha/2);    
else % Exact, binomial law based CI
    CI_lo = max(0, binoinv(alpha/2,n,powertail)/n);
    CI_hi = binoinv(1 - alpha/2,n,powertail)/n;
end

% Plot the results
if plotflag
  loglog(x,1-y,'.')
  hold on
  loglog(x(end-k:end),powertail,'r')
  loglog(x(end-k:end),CI_lo,'r--')
  legend('Data', [num2str(coeff,2) 'x^{' num2str(power,2) '}'], [num2str((1-alpha)*100,2) '% CI'])
  loglog(x(end-k:end),CI_hi,'r--')
  hold off
end

function [POWER,COEFF]=powerfit(x,y)
%POWERFIT Fit power law data.
%   [POWER,COEFF]=POWERFIT(X,Y) finds the coefficients of a power law 
%   Y=COEFF*X^POWER that fits the data (X,Y) in a least squares sense.

x = x(y~=0);
y = y(y~=0);
p = polyfit(log(x(:)),log(y(:)),1);
COEFF = exp(p(2));
POWER = p(1);

