% repetition_curves

% Type:         Script
% Level:        1
% Group:        A
% Usage:        Generates repetition curves
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupA','var')
    error("Run 'A_0_process_data_modelfree.mat' first!");
end

nrnds       = 5;
nrnds_tot   = 2*nrnds; % total number of rounds per condition
nt_pre      = 2;
nt_post     = 16;
repe_cond   = nan(nrnds_tot,nt_pre+nt_post-1,3,nsubj); % total number of rounds per condition, 18 trials, 3 conditions, nsubj

isubj = 0;
for isubj_raw = 1:nsubj_raw
    if ismember(isubj_raw,exclusion_list)
        continue
    else
        isubj = isubj+1;
    end
    
    % identify switch points 
    ind_switch_pts = [0 idx_rond(isubj,2:end)-idx_rond(isubj,1:end-1)];
    idx_switch_pts = find(ind_switch_pts==1);
    idx_cond_start = [1 find(ind_switch_pts==-4)];
    
    % 1/ fill in the 1st rounds for all conditions
    ctr_round = [1 1 1];
    for i = 1:numel(idx_cond_start)
        len_round = min(idx_switch_pts(1+4*(i-1))-idx_cond_start(i),nt_post);
        
        trialrange  = idx_cond_start(i)+1:idx_cond_start(i)+len_round-1; % 1-nt_post
        irond       = ctr_round(idx_cond(isubj,idx_cond_start(i)));
        icond       = idx_cond(isubj,idx_cond_start(i));
        
        % should compare the repetition of [2:nt_post] to [1:nt_post-1] which is then nt_post-1 trials
        repe_cond(irond,3:nt_pre+len_round-1,icond,isubj) = idx_symb(isubj,trialrange) == idx_symb(isubj,trialrange-1);
        ctr_round(icond) = ctr_round(icond)+1;
    end
    
    % 2/ fill in the 0 to 16 trials for REF and UNP and -2 to 16 trials for VOL
    for i = 1:numel(idx_switch_pts)
        % treat the last switch point of each block (indices are multiples of 4)
        if mod(i,4) == 0 
            if i ~= numel(idx_switch_pts)
                % special case for the final switch point of the task
                len_round = min(idx_cond_start(i/4+1)-idx_switch_pts(i),nt_post);
            else
                % last switch point of a block
                len_round = min(size(idx_blck,2)-idx_switch_pts(i)+1,nt_post);
            end
        else
            % intermediate switch points
            len_round = min(idx_switch_pts(i+1)-idx_switch_pts(i),nt_post);
        end
        
        trialrange  = idx_switch_pts(i)-2+1:idx_switch_pts(i)+len_round-1;
        irond       = ctr_round(idx_cond(isubj,idx_switch_pts(i)));
        icond       = idx_cond(isubj,idx_switch_pts(i));
        if icond == 2
            repe_cond(irond,1:nt_pre+len_round-1,icond,isubj) = idx_symb(isubj,trialrange) == idx_symb(isubj,trialrange-1);
        else
            repe_cond(irond,3:nt_pre+len_round-1,icond,isubj) = idx_symb(isubj,trialrange(3:end)) == idx_symb(isubj,trialrange(2:end-1));
        end
        ctr_round(icond) = ctr_round(icond)+1;
    end
end
prepe_cond = permute(squeeze(mean(repe_cond,1,'omitnan')),[2 1 3]);

%plot
colorrgb = [41 52 29; 19 46 55; 58 22 22]/100;
figure;
title(sprintf('Average p(repeat) after round start (all blocks)\nShaded area SEM\nnsubj = %d',nsubj),'FontSize',24);
hold on
for icond = 1:3
    shadedErrorBar(1:nt_pre+nt_post-1,mean(prepe_cond(icond,:,:),3,'omitnan'),std(prepe_cond(icond,:,:),0,3,'omitnan')/sqrt(nsubj),...
                   'lineprops',{'LineWidth',2,'Color',colorrgb(icond,:)});
end
legend({'REF','VOL','UNP'},'Location','southeast')
xticks([1 2 3 5 8 13 17]);
xline(1.5,'HandleVisibility','off');
yline(.5,':','HandleVisibility','off');
xticklabels([-1 1 2 4 8 12 16]);
xlabel('Trial position around round start','FontSize',12);
ylabel('p(repeat)','FontSize',12);
set(gca,'FontSize',12);
set(gca,'TickDir','out');
hold off
