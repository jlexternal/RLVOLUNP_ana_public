% process_data_modelfree

% Type:         Script
% Level:        0 
% Group:        A
% Usage:        This script must be executed before any of the scripts in Group 1 can be run. 
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

clear all
load('../participant_data.mat'); % subj_struct

nsubj_raw = size(subj_struct,2);
ntrls = 80;
nblks = 6;

idx_blck = nan(nsubj_raw,480); % block indices
idx_corr = nan(nsubj_raw,480); % correct reponses indicator
idx_rond = nan(nsubj_raw,480); % round indices
idx_cond = nan(nsubj_raw,480); % condition indices
idx_symb = nan(nsubj_raw,480); % symbol indices
idx_trl  = nan(1,480); % trial number
fb_seen  = nan(nsubj_raw,480);

exclusion_list = [];

% localize data from the table structure
for isubj = 1:nsubj_raw
    resptable = subj_struct{isubj};
    
    idx_blck(isubj,:) = table2array(resptable(:,'i_block'));
    idx_corr(isubj,:) = table2array(resptable(:,'is_correct'));
    idx_rond(isubj,:) = table2array(resptable(:,'i_round'));
    idx_symb(isubj,:) = table2array(resptable(:,'choice_symbol'));
    idx_trl(1,:)      = table2array(resptable(:,'i_trial'));
    fb_seen(isubj,:)  = table2array(resptable(:,'seen_feedback'));

    cond_array = cell2mat(table2array(resptable(:,'cond')));
    cond_array = cond_array(:,1);
    idx_cond(isubj,cond_array=='R') = 1;
    idx_cond(isubj,cond_array=='V') = 2;
    idx_cond(isubj,cond_array=='U') = 3;
    
    if sum(idx_corr(isubj,:)) < 258
        exclusion_list(end+1) = isubj;
    end
end
excl = exclusion_list;
save('../out/excl.mat','excl');
fprintf('File excl saved to /out !\n');
% number of subjects after exclusion criteria
nsubj = nsubj_raw - numel(exclusion_list);

clearvars cond_array resptable

% distribution of feedback value in the course of a round
idx_tr_rel = nan(nsubj,480); % relative trial number (within round)

trl_lens = 1:1:8;

n_lowfb = nan(nsubj,numel(trl_lens)+1,2); % subj, trial length, cond

isubj = 0;
for isubj_raw = 1:nsubj_raw
    if ismember(isubj_raw,exclusion_list)
        continue
    else
        isubj = isubj+1;
    end
    
    ind_samerond = [0 idx_rond(isubj_raw,2:end)-idx_rond(isubj_raw,1:end-1)];
    ind_samerond = ind_samerond == 0;
    ind_samerond(idx_cond(isubj_raw,2:size(idx_cond,2))==2) = 1; % all VOL condition trials are considered the same round
    
    % find index where condition changes
    idx_changecond = [0 idx_cond(isubj_raw,2:end)-idx_cond(isubj_raw,1:end-1)];
    idx_changecond = idx_changecond ~= 0;
    idx_changecond = find(idx_changecond == 1);
    
    % find indices where the round changes (REF and UNP only)
    idx_changerond = find(ind_samerond == 0);
    % find union of the change points
    idx_changepts = union(idx_changecond,idx_changerond);
    
    % relativise trial numbers by round
    for i = 1:numel(idx_changepts)
        idx = idx_changepts(i);
        subtrahend = (idx-1) - floor(idx/80)*80;
        if i == 1
            idx_tr_rel(isubj,1:idx-1) = idx_trl(1:idx-1);
        end
        if i < numel(idx_changepts)
            idx_tr_rel(isubj,idx:idx_changepts(i+1)-1) = idx_trl(idx:idx_changepts(i+1)-1)-subtrahend;
        else
            idx_tr_rel(isubj,idx:end) = idx_trl(idx:end)-subtrahend;
        end
    end
    
    fb = fb_seen(isubj_raw,:);
    for icond = [1 3]
        for i = 1:numel(trl_lens)
            trl_len = trl_lens(i);
            n_lowfb(isubj,i,find([1 3]==icond)) = sum(fb<50 & idx_tr_rel(isubj,:)<=trl_len & idx_cond(isubj_raw,:) == icond);
        end
        
         n_lowfb(isubj,end,find([1 3]==icond)) = sum(fb<50 & idx_cond(isubj_raw,:) == icond);
    end
    
end
clearvars fb i icond ind_samerond idx idx_changecond idx_changepts idx_changerond idx_trl isubj isubj_raw subtrahend subj_struct trl_len trl_lens
isGroupA = true;
fprintf('Processing model-free data (Group A) finished.\n')