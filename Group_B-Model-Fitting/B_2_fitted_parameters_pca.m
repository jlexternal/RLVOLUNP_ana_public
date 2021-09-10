% fitted_parameters_pca

% Type:         Script
% Level:        2
% Group:        B
% Usage:        Run PCA on fitted parameters
% Study:        RLVOLUNP
%
% Outputs:      1/ Matrix of coeffients of the principal components
%               2/ Coefficient values of ingredients (parameters) for each principal
%                   component
%               3/ Scree plot of the principle components and their variance
%                   explained
%               4/ Plot of parameter values split on median PC values
%               5/ Plot of variance explained on INDIVIDUAL conditions from PCs from ALL conditions
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('isGroupB','var')
    error("Run 'B_1_process_model_fits.m' first!");
end

% In running the PCA, we want participants to have understood all aspects of the
%   task. Therefore, we exclude those participants who failed to achieve above chance
%   performance in any condition.
i_excl = [];
i_excl = [i_excl excl];
for icond = 1:3
    i_excl = [i_excl excl_by_cond{icond}'];
end
i = setdiff(idx_subj,i_excl);

ncond         = size(out_pars,3);
npar          = size(out_pars,2);
npar_allpca   = 12;
parstrpca     = {'alpha(R)','delta(R)','zeta(R)','tau(R)','alpha(V)','delta(V)','zeta(V)','tau(V)','alpha(U)','delta(U)','zeta(U)','tau(U)'};
pars_fullexcl = [out_pars(i,:,1) out_pars(i,:,2) out_pars(i,:,3)];

% pca on data
[cpca,spca,lpca,~,epca] = pca(zscore(pars_fullexcl,[],1)); % nsubj, 12 (parameters)
cpca_data = cpca;
spca_data = spca;
save('../out/out_pca_test.mat','spca','cpca');

% sign by 1st ingredient as positive
ind_negalpha = cpca(1,:) < 0;
ind_negalpha = repmat(ind_negalpha,[npar_allpca 1]);
cpca(ind_negalpha) = -cpca(ind_negalpha);

% plot pc coefficient matrix
figure
imagesc(cpca);
colormap('parula');
xticks(1:npar_allpca);
yticks(1:npar_allpca);
xtickstr = cell(npar_allpca,1);
for ipar = 1:npar_allpca
    xtickstr{ipar} = sprintf('%d (%.02f%%)',ipar,epca(ipar));
    for jpar = 1:npar_allpca
        text(ipar,jpar,num2str(cpca(jpar,ipar)),'HorizontalAlignment','center','FontSize',14);
    end
end
colorbar;
xlabel('principle component i (variance explained)');
set(gca,'CLim',[-.5 .5])
set(gca,'XTickLabels',xtickstr);
set(gca,'YTickLabels',parstrpca);
set(gca,'FontSize',10);
title('PCA coefficient matrix');

% bootstrap pca: pick the bootstrapped pc that correlates strongest w/ the original ith PC, and store it
nb  = 1000;
npc = 5;
coeff_bs = nan(npar_allpca,2,nb); % parameter, pcs, bootstrap samples
expl_bs  = nan(npc,nb); % pcs, bootstrap samples
ns       = numel(i);
idx_subj = 1:ns;

% bootstrap samples
fprintf('Bootstrap PCA...\n');
for ipc = 1:npc
    for ib = 1:nb
        idx_bs = randsample(idx_subj,ns,true);
        [coeff,~,~,~,expl] = pca(zscore(pars_fullexcl(idx_bs,:),[],1));
        rhos = corr(cpca(:,ipc),coeff); % correlate w/ ith PC from data
        [~,imax] = max(abs(rhos));
        if rhos(imax) < 0
            coeff(:,imax) = -coeff(:,imax);
        end
        coeff_bs(:,ipc,ib) = coeff(:,imax);
        expl_bs(ipc,ib)    = expl(imax);
    end
end
ind_negalpha = coeff_bs(1,:,:) < 0;
ind_negalpha = repmat(ind_negalpha,[npar_allpca 1 1]);
coeff_bs(ind_negalpha) = -coeff_bs(ind_negalpha);

% 2/ plot coefficients
colorrgb = [41 52 29; 19 46 55; 58 22 22]/100;
figure
for ipc = 1:2
    subplot(2,1,ipc)
    title(sprintf('Value of parameter coefficient (PC%d)\nalpha(R) fixed positive',ipc));
    hold on
    
    for icond = 1:ncond
        idx_cond = [1:4]+(icond-1)*4;
        scatter(idx_cond,cpca(idx_cond,ipc),70,colorrgb(icond,:),'filled');
        scatter(idx_cond,mean(coeff_bs(idx_cond,ipc,:),3),100,colorrgb(icond,:),'x');
        errorbar(idx_cond,mean(coeff_bs(idx_cond,ipc,:),3),std(coeff_bs(idx_cond,ipc,:),0,3),'LineWidth',2,'CapSize',0,'LineStyle','none','Color',colorrgb(icond,:));
        xline(4.5+(icond-1)*4,'LineWidth',1.5);
    end
    
    ylabel('Coefficient value');
    yline(0,':','LineWidth',2);
    xticks(1:npar_allpca);
    xticklabels(parstrpca);
    xlim([.5 npar_allpca+.5]);
end

% validate up to the ith principle component w/ shuffled matrices in pca to find out whether the 
%   var explained is by chance or not (n=1000)
nb = 5000;
coeff_shuff = nan(npar_allpca,npar_allpca,nb); % parameter, pc, bootstrap sample
expl_shuff  = nan(npar_allpca,nb); % pc, bootstrap sample

% shuffle parameter arrays
for ib = 1:nb 
    pars_shuffled = nan(size(pars_fullexcl));
    for ipar = 1:npar_allpca
        idx_bs = randsample(idx_subj,ns,false);
        pars_shuffled(:,ipar) = pars_fullexcl(idx_bs,ipar);
    end
    [coeff_shuff(:,:,ib),~,~,~,expl_shuff(:,ib)] = pca(zscore(pars_shuffled,[],1));
end

% 3/ plot variance explained
figure
sgtitle(sprintf('Percent variance explained\np denotes probability of chance explanation'))
hold on
ymax = 0;
for ipc = 1:npc
    p = mean(expl_shuff(ipc,:) >= epca(ipc));
    % data pcs
    bar(ipc,epca(ipc)/100);
    ymax = max(ymax,epca(ipc)/100+.01);
    text(ipc,epca(ipc)/100+.01,sprintf(' p=%.03f',p),'FontSize',12);
    % shuffled pcs
    errorbar(ipc,mean(expl_shuff(ipc,:))/100,std(expl_shuff(ipc,:)/100,0),'xk','LineWidth',1.5,'CapSize',0);
    % bootstrapped pcs
    scatter(ipc,mean(expl_bs(ipc,:))/100,50,'filled','k')
    errorbar(ipc,mean(expl_bs(ipc,:))/100,std(expl_shuff(ipc,:)/100,0),'k','LineWidth',2,'CapSize',0);
end
ylim([0 ymax+.1]);
xticks(1:npc);
xlabel('Principle component i');
ylabel('Fraction of variance explained');

% 4/ fit parameters split on median PC value
colorrgb    = [41 52 29; 19 46 55; 58 22 22]/100;
nsubj       = size(out_pars,1);
for icond = 1:3
    for ipar = 1:npar
        bar_dat(ipar,icond) = nanmean(out_pars(:,ipar,icond));
        err_dat(ipar,icond) = nanstd(out_pars(:,ipar,icond));
    end
end

figure
sgtitle('Parameter values split on median PC values')
for ipc = 1:2
    for ipar = 1:npar
        subplot(2,4,(ipc-1)*4+ipar);
        title(sprintf('%s \n(median split on PC%d)',parstr{ipar},ipc))
        hold on
        % plot overall mean
        b = bar(bar_dat(ipar,:));
        b.EdgeColor = 'none';
        b.FaceColor = 'flat';
        b.FaceAlpha = .2;
        for icond = 1:3
            b.CData(icond,:) = colorrgb(icond,:);
        end
        % plot means split on PCi
        for icond = 1:3
            e = errorbar(icond,bar_dat(ipar,icond),err_dat(ipar,icond)/sqrt(nsubj),...
                    'Color',colorrgb(icond,:),'LineWidth',1.5,'LineStyle','none','CapSize',0,'HandleVisibility','off');
            set(e.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*.2])
            i = spca(:,ipc) <= nanmedian(spca(:,ipc));
            scatter(icond, mean(out_pars(i,ipar,icond)),50,'v','MarkerFaceColor',colorrgb(icond,:),'MarkerEdgeColor','none');
            errorbar(icond,mean(out_pars(i,ipar,icond)),std(out_pars(i,ipar,icond))/sqrt(nnz(i)),...
                    'Color',colorrgb(icond,:),'LineWidth',1.5,'LineStyle','none','CapSize',0,'HandleVisibility','off');
            i = spca(:,ipc) > nanmedian(spca(:,ipc));
            scatter(icond, mean(out_pars(i,ipar,icond)),50,'^','MarkerFaceColor',colorrgb(icond,:),'MarkerEdgeColor','none');
            errorbar(icond,mean(out_pars(i,ipar,icond)),std(out_pars(i,ipar,icond))/sqrt(nnz(i)),...
                    'Color',colorrgb(icond,:),'LineWidth',1.5,'LineStyle','none','CapSize',0,'HandleVisibility','off');
        end
        xticks(1:3)
        xticklabels(condstr)
        if ipar == 1
            ylabel('Parameter value')
        end
        set(gca,'FontSize',12);
    end
end


% 5/ variance explained on INDIVIDUAL conditions from PCs from ALL conditions

% standardize parameter values
xpar_data = zscore(pars_fullexcl,[],1);

% bootstrap pca: pick the bootstrapped pc that correlates strongest w/ the original ith PC, and store it
nb  = 1000;
npc = 5;
coeff_bs    = nan(npar_allpca,3,nb); % parameter, pcs, bootstrap samples
vexp_cnd_bs = nan(nb,3,npc); 
ns       = numel(i);
idx_subj = 1:ns;
% bootstrap samples
fprintf('Bootstrap PCA...\n');
for ipc = 1:npc
    for ib = 1:nb
        idx_bs  = randsample(idx_subj,ns,true);             % bootstrap sample subject indices
        xpar    = zscore(pars_fullexcl(idx_bs,:),[],1);     % zscore the parameters
        [cpca,spca,lpca,epca] = pca(xpar);                  % run PCA
        rhos = corr(cpca_data(:,ipc),cpca);                 % correlate w/ ith PC from data
        [~,imax] = max(abs(rhos));
        if rhos(imax) < 0
            coeff(:,imax) = -coeff(:,imax);
        end
        coeff_bs(:,ipc,ib) = coeff(:,imax);
        xhat = spca(:,imax)*cpca(:,imax)';
        % analyze the variables for a single component
        xhat    = spca(:,imax)*cpca(:,imax)';
        xres    = xhat-xpar;
        % compute the variance explained for each condition
        vexp_cnd_bs(ib,:,ipc) = [ ...
            1-sum(var(xres(:,1:4)))/sum(var(xpar(:,1:4))), ...
            1-sum(var(xres(:,5:8)))/sum(var(xpar(:,5:8))), ...
            1-sum(var(xres(:,9:12)))/sum(var(xpar(:,9:12)))];
    end
end

% shuffle PCA inputs for stats
vexp_cnd_sh = nan(nb,3,npc);
fprintf('Shuffling PCA inputs...\n');
for ib = 1:nb
    pars_shuffled = nan(size(pars_fullexcl));
    for ipc = 1:npar_allpca
        idx_sh               = randsample(idx_subj,ns,false);
        pars_shuffled(:,ipc) = pars_fullexcl(idx_sh,ipc);
    end
    xpar             = zscore(pars_shuffled,[],1);     % zscore the parameters
    [cpca,spca,lpca] = pca(xpar);                       % run PCA
    for ipc = 1:3
        xhat = spca(:,ipc)*cpca(:,ipc)';
        % analyze the variables for a single component
        xhat    = spca(:,ipc)*cpca(:,ipc)';
        xres    = xhat-xpar;
        % compute the variance explained for each condition
        vexp_cnd_sh(ib,:,ipc) = [ ...
            1-sum(var(xres(:,1:4)))/sum(var(xpar(:,1:4))), ...
            1-sum(var(xres(:,5:8)))/sum(var(xpar(:,5:8))), ...
            1-sum(var(xres(:,9:12)))/sum(var(xpar(:,9:12)))];
    end
end

% plot variance explained per condition
for ipc = 1:2
    xhat_data    = spca_data(:,ipc)*cpca_data(:,ipc)';
    xres_data    = xhat_data-xpar_data;
    vexp_all = 1-sum(var(xres_data))/sum(var(xpar_data));
    % compute the variance explained for each condition
    vexp_cnd = [ ...
        1-sum(var(xres_data(:,1:4)))/sum(var(xpar_data(:,1:4))), ...
        1-sum(var(xres_data(:,5:8)))/sum(var(xpar_data(:,5:8))), ...
        1-sum(var(xres_data(:,9:12)))/sum(var(xpar_data(:,9:12)))];
    figure
    title(sprintf('variance explained on INDIVIDUAL conditions \nby PC%d (from all conditions)',ipc));
    hold on
    b = bar(vexp_cnd);
    for icond = 1:3
        ymax = 0;
        b.FaceColor      = 'flat';
        b.CData(icond,:) = colorrgb(icond,:);
        p = mean(vexp_cnd_sh(:,icond,ipc) >= vexp_cnd_bs(:,icond,ipc));
        % bootstrap shuffled
        scatter(icond,mean(vexp_cnd_sh(:,icond,ipc),1),100,colorrgb(icond,:),'x','filled','MarkerEdgeColor',[1 1 1],'LineWidth',2);
        errorbar(icond,mean(vexp_cnd_sh(:,icond,ipc),1),std(vexp_cnd_sh(:,icond,ipc),1,1),'LineWidth',1,...
                    'CapSize',0,'Color',[1 1 1]);
        % bootstrap
        scatter(icond,mean(vexp_cnd_bs(:,icond,ipc),1),100,colorrgb(icond,:),'filled','MarkerEdgeColor',[1 1 1],'LineWidth',2);
        e = errorbar(icond,mean(vexp_cnd_bs(:,icond,ipc),1),std(vexp_cnd_bs(:,icond,ipc),1,1),'k','LineStyle','none','LineWidth',2,...
                    'CapSize',0,'Color',colorrgb(icond,:).*.7);
        ymax = max(ymax,e.YData+e.YPositiveDelta);
        text(icond,ymax,sprintf(' p=%.03f',p),'FontSize',12);
    end
    xlabel('Condition');
    ylabel('Var explained');
    xticks(1:3);
    xticklabels(condstr);
    set(gca,'FontSize',12);
    ylim([0,0.35]);
end