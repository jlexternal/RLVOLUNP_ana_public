% overall_prepeat_plot

% Type:         Script
% Level:        1
% Group:        A
% Usage:        Generates overall repetition tendency plots
% Study:        RLVOLUNP
%
% Jun Seok Lee  <jlexternal@gmail.com> - September 2021

if ~exist('isGroupA','var')
    error("Run 'A_0_process_data_modelfree.mat' first!");
end

ncond       = 3;
nrepeat     = zeros(nsubj,ncond);
ntotal      = zeros(nsubj,ncond);
nrepeat_raw = zeros(nsubj_raw,ncond);
ntotal_raw  = zeros(nsubj_raw,ncond);

isubj = 0;
for isubj_raw = 1:nsubj_raw
    if ismember(isubj_raw,exclusion_list)
        continue
    else
        isubj = isubj+1;
    end
    
    for iblck = 1:nblks
        icond      = idx_cond(isubj_raw,1+80*(iblck-1));
        ind_blck   = idx_blck(isubj_raw,:) == iblck; 
        symb_blck  = idx_symb(isubj_raw,ind_blck);
        ind_repeat = symb_blck(2:end) == symb_blck(1:end-1);
        nrepeat(isubj,icond) = nrepeat(isubj,icond) + sum(ind_repeat);
        nrepeat_raw(isubj_raw,icond) = nrepeat_raw(isubj_raw,icond) + sum(ind_repeat);
        
        if icond == 2
            ntotal(isubj,icond) = ntotal(isubj,icond) + numel(ind_repeat);
            ntotal_raw(isubj_raw,icond) = ntotal_raw(isubj_raw,icond) + numel(ind_repeat);
        else
            ntotal(isubj,icond) = ntotal(isubj,icond) + numel(ind_repeat) - 4; % -4 for the start of the 4 rounds remaining after the 1st
            ntotal_raw(isubj_raw,icond) = ntotal_raw(isubj_raw,icond) + numel(ind_repeat) - 4;
        end
    end
end

prepeat  = nrepeat./ntotal;
prep_raw = nrepeat_raw./ntotal_raw;

% stats
[p,anova] = anova1(prepeat,[],'off');
if anova{2,6} < .001
    pstr = '<.001';
else
    pstr = sprintf('=%.03f',anova{2,6});
end
anovastr = sprintf('F(%d, %d) = %.02f, p%s',anova{2,3},anova{3,3},anova{2,5},pstr);

% plot
colorrgb = [41 52 29; 19 46 55; 58 22 22]/100;
figure
title(sprintf('Box plot of mean p(repeat)\nANOVA: %s\nDots/lines: Mean/SEM\nPost-hoc tests are sign-rank\nN=%d',anovastr,nsubj));
hold on
bp = boxplot(prepeat,'Notch','on','Labels',{'REF','VOL','UNP'},'Colors',colorrgb);
set(bp,{'linew'},{2})
set(gca,'FontSize',12)
errorbar(1:3,mean(prepeat,1),std(prepeat,0,1)/sqrt(nsubj),'k','LineWidth',2,'LineStyle','none','CapSize',0);
scatter(1:3,mean(prepeat,1),50,colorrgb,'filled');

% post-hoc paired t=tests
cond_pairs = nchoosek(1:3,2);
cond_pairs = [cond_pairs(1,:);cond_pairs(3,:);cond_pairs(2,:)];
for ipair = 1:size(cond_pairs,1)
    pair = cond_pairs(ipair,:);
    x = prepeat(:,pair(1));
    y = prepeat(:,pair(2));
    [~,p] = ttest(x,y);
    if p < .001
        ss = '***';
    elseif p < .01
        ss = '**';
    elseif p < .05
        ss = '*';
    else
        ss = '';
    end
    line([pair(1) pair(2)],.02*ipair+1*ones(2,1),'Color','k');
    text(mean(pair),.02*ipair+1.01,ss,'FontSize',14);
end
ylabel('p(repeat)');
ylim([.37 1.1])
hold off;
