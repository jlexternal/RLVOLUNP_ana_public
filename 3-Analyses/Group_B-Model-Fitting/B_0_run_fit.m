% run_fits

% Type:         ??
% Level:        0
% Group:        B
% Usage:        Fit the model to data
% Study:        RLVOLUNP
%
% ???

% fit_noisyKF_cfrule.m is the fitting code
%   our settings:
%       cfrule  : true
%       nstype  : 'weber'
%       chrule  : 'softm'
%       nsmp    : 1e3

% Whatever optimizer you use to fit the model, the fitted parameter outputs must
%   be organized in the matrix named 'pars' with the following dimensions :
%   (n participants x n parameters x n condition)
% 
% For example:
%   pars = nan(200,4,3);


% --- Write your model fitting script here. ---