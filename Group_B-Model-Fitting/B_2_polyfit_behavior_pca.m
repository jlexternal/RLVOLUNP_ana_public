% polyfit_behavior_pca

% Type:         Script
% Level:        2
% Group:        B
% Usage:        Fit polynomial model to performance vs PC score data
% Study:        RLVOLUNP
%
% Requirements: This script requires A_1_overall_accuracy_plot to be run first for
%                   the pcor_raw.mat file to be stored in the /out/ folder.
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupB','var')
    error("Run 'B_1_process_model_fits.m' first!");
end
addpath('../functions/');
addpath('../out/');

if ~exist('pcor_raw.mat','file')
    error("Run 'A_1_overall_accuracy_plot.m' first to generate pcor_raw.mat file!");
end
load('pcor_raw.mat'); % loads pcor_raw

i_excl = [];
i_excl = [i_excl excl];
for icond = 1:3
    i_excl = [i_excl excl_by_cond{icond}'];
end
i = setdiff(idx_subj,i_excl);
    
% run PCA
allpars  = [out_pars(i,:,1) out_pars(i,:,2) out_pars(i,:,3)];
[cpca,spca,lpca,~,epca] = pca(zscore(allpars,[],1));

% fit p(cor) as quadratic function of PCi
npc         = 2;
ncond       = size(out_pars,3);
condstr     = {'REF','VOL','UNP'};
modstr      = {'linear','quadratic'};
colorrgb    = [41 52 29; 19 46 55; 58 22 22]/100;
r_cor       = cell(npc,ncond,3);    % pcs, conds, polydeg
f_cor       = cell(npc,ncond,3);    % pcs, conds, polydeg
b_cor       = nan(4,npc,ncond,3);   % coefs, pcs, conds, polydeg
pv_cor      = nan(4,npc,ncond,3);   % pvals, pcs, conds, polydeg
bic_cor     = nan(npc,ncond,3);     % pcs, conds, polydeg

% fit NLMs
for ipc = 1:npc
    for icond = 1:ncond
        X = spca(:,ipc);
        y = pcor_raw(i,icond);
        for ideg = 1:2 % degrees of the polynomial
            f = @(b,x) polyval_asc(b,x);    % fitting model
            if ideg == 1
                r = fitlm(X,y);   
                f = fit(X,y,'poly1');
            else
                r = fitnlm(X,y,f,zeros(1,ideg+1));   % result
                f = fit(X,y,'poly2');
            end
            r_cor{ipc,icond,ideg}           = r;
            f_cor{ipc,icond,ideg}           = f;
            b_cor(1:ideg+1,ipc,icond,ideg)  = r.Coefficients{:,1};
            pv_cor(1:ideg+1,ipc,icond,ideg) = r.Coefficients{:,4};
            bic_cor(ipc,icond,ideg)         = r.ModelCriterion.BIC;
        end
    end
end

% plot winning models
coeffs = nan(3,3,npc); % condition, coefficient, ith pc
figure
sgtitle(sprintf('Fit of p(cor) ~ Pn(pc_i)','Interpreter','latex'))
for ipc = 1:npc
    for icond = 1:ncond
        [~,imodel] = min(squeeze(bic_cor(ipc,icond,:)));
        subplot(2,3,(ipc-1)*3+icond)
        hold on
        title(sprintf('Condition: %s\n model: %s',condstr{icond},modstr{imodel}))
        X = spca(:,ipc);
        y = pcor_raw(i,icond);
        r = r_cor{ipc,icond,imodel};
        
        scatter(X,y,40,colorrgb(icond,:),'filled');
        xrange  = (min(X)*1.1):.1:(max(X)*1.1);
        pn      = polyval_asc(r.Coefficients{:,1}',xrange);
        if ipc == 1 
            coeffs(icond,1:2,ipc) = r.Coefficients{:,1}';
        else
            coeffs(icond,:,ipc) = r.Coefficients{:,1}';
        end
        b_er    = coefCI(r);
        if size(b_er,2) == 1
            pn_ci_up  = polyval_asc(b_er(:,1),xrange);
            pn_ci_dn  = polyval_asc(b_er(:,2),xrange);
        else
            pn_ci_up  = polyval_asc(b_er(:,1)',xrange);
            pn_ci_dn  = polyval_asc(b_er(:,2)',xrange);
        end
        plot(xrange,pn)
        s = shadedErrorBar(xrange,pn,[abs(pn-pn_ci_up); abs(pn-pn_ci_dn)],'patchSaturation',.1,...
                        'lineprops',{'LineWidth',2,'Color',colorrgb(icond,:)});
        set(s.edge,'LineStyle','none');
        
        ep = predint(f_cor{ipc,icond,imodel},xrange,.95);
        sp = shadedErrorBar(xrange,pn,[abs(pn-ep(:,2)'); abs(pn-ep(:,1)')],'patchSaturation',.1,...
                    'lineprops',{'LineWidth',.5,'Color',colorrgb(icond,:)});
        set(sp.edge,'LineWidth',2,'LineStyle',':');
        
        if icond == 1
            ylabel('p(cor)');
        elseif icond == 2
            xlabel(sprintf('z-score (PC%d)',ipc));
        end
        set(gca,'FontSize',12);
    end
end

% plot areas of confidence interval (or lack of) overlap
figure
sgtitle(sprintf('Confidence intervals of the fitted polynomial models\n95%% range of data shown'));
hold on
for ipc = 1:npc
    subplot(2,1,ipc)
    X_agg = [];
    X_avg = nan(3,1);
    for icond = 1:ncond
        [~,imodel] = min(squeeze(bic_cor(ipc,icond,:)));
        hold on
        title(sprintf('PC%d\nWinning model: %s',ipc,modstr{imodel}))
        X               = spca(:,ipc);
        X_agg           = [X_agg; X];
        X_avg(icond)    = mean(X);
        y               = pcor_raw(i,icond);
        r               = r_cor{ipc,icond,imodel};
        
        xrange  = (min(X)*1.1):.1:(max(X)*1.1);
        pn      = polyval_asc(r.Coefficients{:,1}',xrange);
        b_er    = coefCI(r);
        if size(b_er,2) == 1
            pn_ci_up  = polyval_asc(b_er(:,1),xrange);
            pn_ci_dn  = polyval_asc(b_er(:,2),xrange);
        else
            pn_ci_up  = polyval_asc(b_er(:,1)',xrange);
            pn_ci_dn  = polyval_asc(b_er(:,2)',xrange);
        end
        
        plot(xrange,pn);
        s = shadedErrorBar(xrange,pn,[abs(pn-pn_ci_up); abs(pn-pn_ci_dn)],'patchSaturation',.4,...
                        'lineprops',{'LineWidth',2,'Color',colorrgb(icond,:)});
        set(s.edge,'LineStyle','none');
        scatter(X_avg(icond),polyval_asc(r.Coefficients{:,1}',X_avg(icond)),40,'MarkerFaceColor',colorrgb(icond,:),'MarkerEdgeColor',[1 1 1],'LineWidth',2);
    end
    xlim([mean(X_agg)-2*std(X_agg) mean(X_agg)+2*std(X_agg)]);
    set(gca,'FontSize',12);
    set(gca,'TickDir','out');
end
xlabel('Z-score(PC)')
hold off

% function that inverts the order of polynomial degree arguments (b)
function y = polyval_asc(b,x)
    y = polyval([b(2:end) b(1)],x);
end