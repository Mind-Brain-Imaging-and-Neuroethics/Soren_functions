function V = vif(data)

R0 = corrcoef(data);
V = diag(inv(R0))';