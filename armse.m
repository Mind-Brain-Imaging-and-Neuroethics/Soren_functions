function [mse] = armse(model,data,varargin)

%model = ar(iddata(vert(data)),order);

mse = model.Report.Fit.MSE;

