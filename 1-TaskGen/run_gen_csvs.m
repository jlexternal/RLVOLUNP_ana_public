% run_gen_csvs

% Type:         Script
% Level:        3
% Group:        TaskGen
% Usage:        Generates CSV files containing task data for a range of participants in the
%                   online task
% Study:        RLVOLUNP
%
% Requirements: betafun.m, betapar.m, betasmp.m, gen_drift.m, gen_sessions.m, gen_csv_for_online_expe.m
%               unique_id_list.mat (from gen_unique_id_list.m)
%
% Usage:        Call function in a simple for loop in the command window iterating over the
%                   subject numbers for whom you wish to generate experiment CSVs.
%
% Note:         Once files with a specific subject number (enumerator_id) has been used,
%                   it should not be used again.
%               
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

% Indicate range:
subjrange = 1:100;

for isubj = subjrange
    fprintf('Generating online task CSV file for participant %d...\n',isubj);
    gen_csv_for_online_expe(isubj);
end