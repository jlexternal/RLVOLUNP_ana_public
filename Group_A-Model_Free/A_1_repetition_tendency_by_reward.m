% repetition_tendency_by_reward

% Type:         Script
% Level:        1
% Group:        A
% Usage:        Generate average repetition tendency by reward magnitude plots
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupA','var')
    error("Run 'A_0_process_data_modelfree.mat' first!");
end

% bin grain
bin_right_lim = 20:20:100;

nrepeat = nan(nsubj,numel(bin_right_lim),2,3); % subject, bin, round split, condition
nswitch = nan(nsubj,numel(bin_right_lim),2,3);

isplit = 5; % inclusive on less than, exclusive on greater than

isubj = 0;
for isubj_raw = 1:nsubj_raw
    if ismember(isubj_raw,exclusion_list)
        continue
    else
        isubj = isubj+1;
    end
    
    % identify repeat trials
    ind_repeat = idx_symb(isubj_raw,2:end)-idx_symb(isubj_raw,1:end-1);
    
    % identify non-repeat trials
    ind_switch = ind_repeat ~= 0;
    ind_repeat = ind_repeat == 0;
    
    % for non-repeat trials, need exclude round-change trials
    ind_samerond = idx_rond(isubj_raw,2:end)-idx_rond(isubj_raw,1:end-1);
    ind_samerond = ind_samerond == 0;
    % (here, we are considering that rounds in VOL are separate)
    
    % identify same condition trials 
    ind_samecond = idx_cond(isubj_raw,2:end)-idx_cond(isubj_raw,1:end-1);
    ind_samecond = ind_samecond == 0;
    
    % find intersection for repeats and switches
    % ___ NOTE: these are trials 2:end ___ %
    ind_repeat = ind_repeat & ind_samecond;
    ind_switch = ind_switch & ind_samecond & ind_samerond;
    
    % identify reward values and bin
    fb = fb_seen(isubj_raw,1:end-1);
    for icond = 1:3
        for ibin = 1:numel(bin_right_lim)
            if ibin == 1
                ind_fb_inrange = fb>0 & fb<=bin_right_lim(1);
            else
                ind_fb_inrange = fb>bin_right_lim(ibin-1) & fb<=bin_right_lim(ibin);
            end
            
            for ihalf = 1:2
                if ihalf == 1
                    idx_split = idx_tr_rel(isubj,2:end) <= isplit;
                else
                    idx_split = idx_tr_rel(isubj,2:end) > isplit;
                end
                % count the number of repeats and switches within feedback bin range
                nrepeat(isubj,ibin,ihalf,icond) = sum(ind_repeat & ind_fb_inrange & idx_cond(isubj_raw,2:end) == icond & idx_split); % responses from feedback within range from chosen condition
                nswitch(isubj,ibin,ihalf,icond) = sum(ind_switch & ind_fb_inrange & idx_cond(isubj_raw,2:end) == icond & idx_split);
            end
        end
    end
end

 % calculate p(repeat) as a function of reward magnitude
prepeat = nrepeat./(nrepeat+nswitch); % subj, bin, split, cond

bin_midpts    = nan(size(bin_right_lim));
bin_midpts(1) = mean([bin_right_lim(1) 0]);
for ib = 2:numel(bin_midpts)
    bin_midpts(ib) = mean([bin_right_lim(ib) bin_right_lim(ib-1)]);
end

% plot
colorrgb    = [41 52 29; 19 46 55; 58 22 22]/100;
halfMarker  = {'x:','o-'};
condstr     = {'REF','VOL','UNP'};

% plots by condition
figure;
hold on
sgtitle(sprintf('p(repeat) as a function of previous reward magnitude split at trial %d\nBins by 20\nnsubj = %d\n(variable SEM)',isplit,nsubj),'FontSize',16);
for icond = 1:3
    subplot(3,1,icond)
    title(condstr{icond})
    
    for ihalf = 1:2
        nsubj_local = sum(~isnan(prepeat(:,:,ihalf,icond)),1);
        shadedErrorBar(bin_midpts,mean(prepeat(:,:,ihalf,icond),1,'omitnan'),std(prepeat(:,:,ihalf,icond),0,1,'omitnan')./sqrt(nsubj_local),...
                       'lineprops',{halfMarker{ihalf},'LineWidth',2,'Color',colorrgb(icond,:)});
    end
    
    legend({sprintf('Trials <=%d',isplit),sprintf('Trials >%d',isplit)},'Location','southeast')
    yline(.5,'--','HandleVisibility','off');
    xlim([0 100])
    
    for binlim = bin_right_lim
        xline(binlim,':','HandleVisibility','off');
    end
    
    xticks(bin_midpts);
    xticklabels({'0-20','21-40','41-60','61-80','81-100'});
    ylabel('p(repeat)');
    xlabel('reward range');
    set(gca,'FontSize',12);
    set(gca,'TickDir','out');
end
hold off

% plots by episode split
figure;
hold on
sgtitle(sprintf('p(repeat) as a function of previous reward magnitude split at trial %d\nBins by 20\nnsubj = %d\n(variable SEM)',isplit,nsubj),'FontSize',16);
for ihalf = 1:2
    subplot(1,2,ihalf)
    
    if ihalf == 1
        title(sprintf('Trials <=%d',isplit));
    else
        title(sprintf('Trials >%d',isplit));
    end
    
    for icond = 1:3
        nsubj_local = sum(~isnan(prepeat(:,:,ihalf,icond)),1);
        shadedErrorBar(bin_midpts,mean(prepeat(:,:,ihalf,icond),1,'omitnan'),std(prepeat(:,:,ihalf,icond),0,1,'omitnan')./sqrt(nsubj_local),...
                       'lineprops',{halfMarker{ihalf},'LineWidth',2,'Color',colorrgb(icond,:)});
    end
    
    legend({'REF','VOL','UNP'},'Location','southeast');
    yline(.5,'--','HandleVisibility','off');
    xlim([0 100])
    
    for binlim = bin_right_lim
        xline(binlim,':','HandleVisibility','off');
    end
    
    xticks(bin_midpts);
    xticklabels({'0-20','21-40','41-60','61-80','81-100'});
    ylabel('p(repeat)');
    xlabel('reward range');
    set(gca,'FontSize',12);
    set(gca,'TickDir','out');
end
hold off
