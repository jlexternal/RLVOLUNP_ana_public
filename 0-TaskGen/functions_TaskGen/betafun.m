function [p] = betafun(x,m,v)
% get pdf of beta distribution with mean m and variance v
x = bsxfun(@plus,x,zeros(size(m)));
[a,b] = betapar(m,v);
p = betapdf(x,a,b);
end