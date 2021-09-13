function [x] = betasmp(m,v)
% sample from beta distribution with mean m and variance v
[a,b] = betapar(m,v);
x = betarnd(a,b);
end