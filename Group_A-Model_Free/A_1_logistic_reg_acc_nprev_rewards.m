% logistic_reg_acc_nprev_rewards

% Type:         Script
% Level:        1
% Group:        A
% Usage:        Regularized logistic regression of accuracy by previous rewards:
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupA','var')
    error("Run 'A_0_process_data_modelfree.mat' first!");
end

addpath ../functions % requires function logreg_regul.m

% logistic regression parameters
bigN    = 5; % max trials back to consider
b_corr  = nan(nsubj,bigN,3);
nrun    = 10; % number of random starting points for regularized logistic regression
s       = 10*ones(1,bigN); % regularizing prior s.d.

isubj = 0;
for isubj_raw = 1:nsubj_raw
    if ismember(isubj_raw,exclusion_list)
        continue
    else
        isubj = isubj+1;
    end
    
    for icond = 1:3
        fb_nback = []; % data structure of rewards 1-n back from choice
        resps    = []; % responses corresponding to set of n-back feedback 
        
        % 1/ go through indices that are less than N and take partial data
        if true
            for i = 2:bigN  % (choices from 1st trial of any round is NOT influenced by any feedback)
                idx = find(idx_tr_rel(isubj,:)==i & idx_cond(isubj_raw,:)==icond); % i trials after start of round
                fb_iback_arr = [];
                
                for j = 1:i-1
                    resp_temp = idx_corr(isubj_raw,idx-j);
                    resp_temp(resp_temp == 0) = -1;
                    
                    fb_temp      = (fb_seen(isubj_raw,idx-j)/100-.5) .* resp_temp;
                    fb_iback_arr = cat(2,fb_iback_arr,fb_temp'); % feedback from idx-j
                end
                
                fb_nback = [fb_nback; cat(2,fb_iback_arr,nan(numel(idx),6-i))];
                resps    = [resps; idx_corr(isubj_raw,idx)']; % resp at idx
            end
        end

        % 2/ identify all indices to take data from (idx > N)
        idx          = find(idx_tr_rel(isubj,:) > bigN & idx_cond(isubj_raw,:) == icond);
        fb_iback_arr = [];
        
        for j = 1:bigN
            resp_temp = idx_corr(isubj_raw,idx-j);
            resp_temp(resp_temp == 0) = -1;
            
            fb_temp      = (fb_seen(isubj_raw,idx-j)/100-.5) .* resp_temp;
            fb_iback_arr = cat(2,fb_iback_arr,fb_temp');
        end
        
        fb_nback     = [fb_nback; fb_iback_arr];
        resps        = [resps; idx_corr(isubj_raw,idx)'];
        
        fprintf('Running logistic regression on participant %d on condition %d...\n',isubj_raw,icond);
        b_corr(isubj,:,icond) = logreg_regul(fb_nback,resps,s,nrun,true);
    end
end

% plot
colorrgb    = [41 52 29; 19 46 55; 58 22 22]/100;

figure
title(sprintf('Regularized logistic regression weights of p(correct) using previous n rewards\nShaded area SEM\nN = %d'...
               ,nsubj),'FontSize',16);
hold on
for icond = 1:3
    shadedErrorBar(1:5,mean(b_corr(:,:,icond)),std(b_corr(:,:,3))/sqrt(nsubj),...
        'lineprops',{'o-','LineWidth',2,'Color',colorrgb(icond,:)});
end
yline(0,'HandleVisibility','off');
legend({'REF','VOL','UNP'},'Location','southwest');
xticks(1:5)
xlabel('n trials back');
xlim([.5 5.5]);
ylabel('logistic weight coefficient');
set(gca,'FontSize',12);
set(gca,'TickDir','out');
hold off