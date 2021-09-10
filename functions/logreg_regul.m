function [b] = logreg_regul(x,y,s,nrun,disc)

if nargin < 3
    s = [];
end
if nargin < 2
    error('Missing input arguments!');
end

% set additional parameters
warnoff = true; % disable warnings?
leps    = 1e-6; % infinitesimal log-evidence

% get mask for missing values
w = ~isnan(x);

% set missing values to zero
x(~w) = 0;

% ensure binary output array
y = logical(y);

% define fitting options
options = optimset('fminsearch');
if warnoff
    options = optimset(options,'Display','off');
else
    options = optimset(options,'Display','notify','MaxFunEvals',1e4);
end

n = size(x,2);
if ~isempty(s) && numel(s) ~= n
    error('Wrong regularizing prior s.d.!');
end
s = s(:);

% fit logistic regression model with regularization
pval = nan(nrun,n);
fval = nan(nrun,1);
for irun = 1:nrun
    % set random starting point
    if ~isempty(s)
        p0 = normrnd(0,s);
    else
        % use uniform prior in [-10,+10]
        p0 = unifrnd(-10,+10,[n,1]);
    end
    [pval(irun,:),fval(irun)] = fminsearch(@(ps)-getl(ps),p0,options);
end

% get best fit
[~,irun] = min(fval);
b = pval(irun,:); % coefficient estimates
l = -fval(irun); % model log-evidence

    function [l] = getl(p)
        % get model log-evidence
        l = normcdf(x*p);
        l(y == 0) = 1-l(y == 0);
        l = sum(log(max(l,leps)));
        if ~isempty(s)
            % add regularizing prior
            if disc % discounted prior
                k = mean(w,1)';
                l = l+sum(k.*log(normpdf(p,0,s)));
            else % fixed prior
                l = l+sum(log(normpdf(p,0,s)));
            end
        end
    end

end