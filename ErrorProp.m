function [delR,R] = ErrorProp(x,delx,y,dely,A,B,R,formType)

%R = form(x,y);

if strcmpi(formType,'add')
    %syms delform(A,delx,B,dely)
    %delform(A,delx,B,dely) = ;
    delR = sqrt((A^2)*(delx^2) + (B^2)*(dely^2));
elseif strcmpi(formType,'multiply')
    %    syms delform(A,delx,B,dely)
    %delform(A,x,delx,B,y,dely,R) = abs(R)*sqrt((A^2)*((delx/x)^2) + (B^2)*((dely/y)^2));
    delR = abs(R)*sqrt((A^2)*((delx/x)^2) + (B^2)*((dely/y)^2));
end
