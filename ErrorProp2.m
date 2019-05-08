function [delR,R] = ErrorProp2(a,dela,b,delb,formula)

syms f(x,y)
syms g(x,y,delx,dely)

f(x,y) = eval(formula);

g(x,y,delx,dely) = sqrt((diff(f,x)^2)*(delx^2) + (diff(f,y)^2)*(dely^2));

R = double(f(a,b));
delR = double(g(a,b,dela,delb));




%R = form(x,y);

% if strcmpi(formType,'add')
%     %syms delform(A,delx,B,dely)
%     %delform(A,delx,B,dely) = ;
%     delR = sqrt((A^2)*(delx^2) + (B^2)*(dely^2));
% elseif strcmpi(formType,'multiply')
%     %    syms delform(A,delx,B,dely)
%     %delform(A,x,delx,B,y,dely,R) = abs(R)*sqrt((A^2)*((delx/x)^2) + (B^2)*((dely/y)^2));
%     delR = abs(R)*sqrt((A^2)*((delx/x)^2) + (B^2)*((dely/y)^2));
% end
