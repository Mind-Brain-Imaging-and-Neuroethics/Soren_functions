function CI = corrci(r,n)

z = 0.5 * log((1 + r) / (1 - r));
se = sqrt((1 / (n - 3)));

CI_z = [z-1.96*se z+1.96*se];

CI = (exp(2*CI_z)-1)./(exp(2*CI_z)+1);


