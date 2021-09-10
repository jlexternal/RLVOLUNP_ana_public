% fit_noisyKF_cfrule

% Type:         Function
% Level:        0
% Group:        B
% Usage:        Samples from truncated normal distribution
% Study:        RLVOLUNP
%

function [x] = noisrnd(m,s)
% sample from noise distribution = truncated normal in [0,1]
if numel(s) == 1
    s = s(ones(size(m)));
end
x = nan(size(m));
f = true(size(m));
while any(f(:))
    x(f) = normrnd(m(f),s(f));
    f(x >= 0 & x <= 1) = false;
end
end