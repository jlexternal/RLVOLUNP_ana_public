% overall_accuracy_plot

% Type:         Script
% Level:    	1
% Group:        A
% Usage:        Generates overall accuracy plots
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupA','var')
    error("Run 'A_0_process_data_modelfree.mat' first!");
end

ncond           = 3;
ncorrect        = zeros(nsubj,ncond);
ntotal          = zeros(nsubj,ncond);
ncorrect_raw    = zeros(nsubj_raw,ncond);
ntotal_raw      = zeros(nsubj_raw,ncond);
condexcl        = zeros(nsubj,ncond);
condexcl_raw    = zeros(nsubj,ncond);
isubj_excl      = cell(3,1);
isubj_raw_excl  = cell(3,1);

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
        corr_blck  = idx_corr(isubj_raw,ind_blck);
        
        ncorrect(isubj,icond)   = ncorrect(isubj,icond) + sum(corr_blck);
        ntotal(isubj,icond)     = ntotal(isubj,icond) + numel(corr_blck); 
        
        ncorrect_raw(isubj_raw,icond) = ncorrect_raw(isubj_raw,icond) + sum(corr_blck);
        ntotal_raw(isubj_raw,icond)   = ntotal_raw(isubj_raw,icond) + numel(corr_blck); 
    end
    
    % id participants that failed below chance on any given condition
    for icond = 1:3
        condexcl(isubj,icond)           = binocdf(ncorrect(isubj,icond),160,.5,'upper') >= .05;
        condexcl_raw(isubj_raw,icond)   = binocdf(ncorrect(isubj,icond),160,.5,'upper') >= .05;
    end
end

% save participants excluded by performance on any one condition to file
for icond = 1:3
    isubj_excl{icond} = find(condexcl(:,icond) == 1);
    isubj_raw_excl{icond} = find(condexcl_raw(:,icond) == 1);
end
excl_by_cond = isubj_raw_excl;
save('../out/excl_by_cond.mat','excl_by_cond');
fprintf('File excl_by_cond saved to /out !\n');

% save average performance
pcorrect = ncorrect./ntotal;
pcor_raw = ncorrect_raw./ntotal_raw;
save('../out/pcor_raw.mat','pcor_raw');
fprintf('File pcor_raw saved to /out !\n');

% stats
[p,anova] = anova1(pcorrect,[],'off');
if anova{2,6} < .001
    pstr = '<.001';
else
    pstr = sprintf('=%.03f',anova{2,6});
end
anovastr = sprintf('F(%d, %d) = %.02f, p%s',anova{2,3},anova{3,3},anova{2,5},pstr);

% plot
colorrgb = [41 52 29; 19 46 55; 58 22 22]/100;
figure
title(sprintf('Box plot of mean p(correct)\nANOVA: %s\nDots/lines: Mean/SEM\nPost-hoc paired t-tests\nN=%d',anovastr,nsubj));
hold on
bp = boxplot(pcorrect,'Notch','on','Labels',{'REF','VOL','UNP'},'Colors',colorrgb);
set(bp,{'linew'},{2})
set(gca,'FontSize',12)
errorbar(1:3,mean(pcorrect,1),std(pcorrect,0,1)/sqrt(nsubj),'k','LineWidth',2,'LineStyle','none','CapSize',0);
scatter(1:3,mean(pcorrect,1),50,colorrgb,'filled');

% post-hoc tests
cond_pairs = nchoosek(1:3,2);
cond_pairs = [cond_pairs(1,:);cond_pairs(3,:);cond_pairs(2,:)];
for ipair = 1:size(cond_pairs,1)
    pair = cond_pairs(ipair,:);
    x = pcorrect(:,pair(1));
    y = pcorrect(:,pair(2));
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
ylabel('p(correct)');
ylim([.37 1.1])
hold off;