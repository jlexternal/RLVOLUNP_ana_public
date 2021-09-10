% process_model_fits

% Type:         Script
% Level:        1
% Group:        B
% Usage:        Process model fit data for Group B analyses
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

% Note: Whatever method you choose to fit the model, the fitted parameter outputs must
%           be organized in the matrix named 'pars' with the following dimensions :
%           (participants x parameter x condition)

% Note: If any participants are to be excluded for any reason, indicate with
%           idx_subj. However, these excluded participants must be accounted for
%           (have a place in the eventual out_pars matrix)

clear all

if ~exist('../out/pars.mat','file')
    error('Matrix containing fitted parameters missing in /out/ directory!');
end

% exclusions for analyses using PCA exclude participants who performed
%   poorly in any one condition
if ~exist('../out/excl.mat','file')
    error('Run A_0_process_data_modelfree.m first to generate overall exclusion list!');
end

if ~exist('../out/excl_by_cond.mat','file')
    error('Run A_1_overall_accuracy_plot.m first to generate exclusion-by-condition list!');
end
addpath('../out/');
load('excl.mat');           % loads excl
load('excl_by_cond.mat');   % loads excl_by_cond

% User modification required: 
% idx_subj = []; % array of participant indices to be included in analysis
idx_subj = 1:200;

parstr  = {'alpha','delta','zeta','tau'};
parrgb  = [80 80 80; 94 69 43; 46 78 93; 81 50 81]/100;
condstr = {'REF','VOL','UNP'};

% temporarily
load ../out/pars.mat
out_pars = pars;

isGroupB = true;