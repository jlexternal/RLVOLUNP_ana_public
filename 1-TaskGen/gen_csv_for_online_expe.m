% gen_csv_for_online_expe

% Type:         Function
% Level:        2
% Group:        TaskGen
% Usage:        Generates CSV files containing task data for one participant in the
%                   online task
% Study:        RLVOLUNP
%
% Requirements: betafun.m, betapar.m, betasmp.m, gen_drift.m, gen_sessions.m
%               unique_id_list.mat (from gen_unique_id_list.m)
%               
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

function gen_csv_for_online_expe(isubj)

iscontinue = true;

if not(isfolder('online_expe_csvs'))
    mkdir online_expe_csvs
end

% check for folder/files already created with with the input isubj
while iscontinue
    if isfolder(sprintf('./online_expe_csvs/subj_%03d',isubj))
        error('Folder with subject number %03d already exists!',isubj);
        iscontinue = false;
    else
        mkdir(sprintf('online_expe_csvs/subj_%03d',isubj));

        % import the unique_id_list
        unique_id_list = load('unique_id_keys/unique_id_list.mat','unique_id_list');
        unique_id_list = unique_id_list.unique_id_list;

        % set global random number stream
        rng = RandStream('mt19937ar','Seed',isubj);
        RandStream.setGlobalStream(rng);

        % generate unique_id_key (key for linking data in the final SQL database)
        unique_id_key = unique_id_list{isubj};

        task_struct = {};
        for i = 1:2
            task_struct{i} = gen_sessions;
        end

        traj_out = [task_struct{1}.traj; task_struct{2}.traj];
        idx_epi = [task_struct{1}.idx_epi; task_struct{2}.idx_epi];

        for i = 1:2
            ind_even = mod(task_struct{i}.idx_epi,2)==0;
            traj_out([1:3]+(3*(i-1)),ind_even==1) = 1-traj_out([1:3]+(3*(i-1)),ind_even==1);
        end
        
        % export data for online experiment
        %   traj_out:   the reward value for the CORRECT OPTION
        %   idx_epi:    block label for each trial (to identify new/switch trials)
        traj_filename       = sprintf('traj_%s.csv',unique_id_key);
        idx_epi_filename    = sprintf('idx_epi_%s.csv',unique_id_key);
        
        csvwrite(traj_filename,traj_out);
        csvwrite(idx_epi_filename,idx_epi);
        
        movefile(traj_filename,sprintf('./online_expe_csvs/subj_%03d/',isubj));
        movefile(idx_epi_filename,sprintf('./online_expe_csvs/subj_%03d/',isubj));
        
        iscontinue = false;
    end
end

end