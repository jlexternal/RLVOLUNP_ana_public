% fitted_parameters_overview

% Type:         Script
% Level:        2
% Group:        B
% Usage:        Run statistics on fitted parameters
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupB','var')
    error("Run 'B_1_process_model_fits.mat' first!");
end

i = setdiff(idx_subj,excl);

% stats on parameters
parvec             = out_pars(i,:,:);
[nsubj,npar,ncond] = size(parvec);

% repeated-measures ANOVA 
fprintf('Running repeated-measures ANOVA on fit parameters...\n');
for ipar = 1:npar
    fprintf('Parameter: %s\n',parstr{ipar});
    repanova_auto(squeeze(parvec(:,ipar,:)),{'condition'});
end

ss = cell(3,npar); % signficance stars (condition combinations, params)
ps = nan(3,npar);  % p-values

% post-hoc ttests
combis = nchoosek(1:ncond,2); % condition pair combinations
for ipar = 1:npar
    for ic = 1:size(combis,1)
        [~,ps(ic,ipar)] = ttest(squeeze(parvec(:,ipar,combis(ic,1))),squeeze(parvec(:,ipar,combis(ic,2))));
        if ps(ic,ipar) < .001
            ss{ic,ipar} = '***';
        elseif ps(ic,ipar) < .01
            ss{ic,ipar} = '**';
        elseif ps(ic,ipar) < .05
            ss{ic,ipar} = '*';
        else
            ss{ic,ipar} = '';
        end
    end
end

% plot
bar_dat = nan(npar,ncond);
err_dat = nan(npar,ncond);
colorrgb    = [41 52 29; 19 46 55; 58 22 22]/100;

figure
clf
hold on
sgtitle(sprintf('Fit parameters\nerror bars SEM\nnsubj=%d',nsubj),'FontSize',16,'Interpreter','none');
for icond = 1:3
    for ipar = 1:npar
        bar_dat(ipar,icond) = mean(parvec(:,ipar,icond),'omitnan');
        err_dat(ipar,icond) = std(parvec(:,ipar,icond),'omitnan');
    end
end

for ipar = 1:npar
    subplot(1,npar,ipar);
    title(parstr{ipar})
    
    hold on
    h = bar(bar_dat(ipar,:));
    % error bars
    for icond = 1:ncond
        errorbar(icond,bar_dat(ipar,icond),err_dat(ipar,icond)/sqrt(nsubj),...
                'Color','k','LineWidth',1.5,'LineStyle','none','CapSize',0,'HandleVisibility','off');
    end
    
    xticks(1:3);
    xticklabels(condstr);
    
    % post-hoc stats
    ymax = max(bar_dat(ipar,:) + err_dat(ipar,:)/sqrt(nsubj),[],2);
    
    for ic = 1:size(combis,1)
        if ~isempty(ss{ic,ipar})
            x = combis(ic,:);
            line(x,ones(2,1)*ymax+.01*(ic),'Color','k');
            text(mean(x),ymax+.011*(ic),ss{ic,ipar},'FontSize',14,'HorizontalAlignment','center');
        end
    end
    hold off
end