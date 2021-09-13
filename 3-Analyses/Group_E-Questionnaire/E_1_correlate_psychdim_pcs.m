% correlate_psychdim_pcs

% Type:         Script
% Level:        1
% Group:        E
% Usage:        Look for correlations between psychiatric dimension scores and
%                   principal components
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

if ~exist('../out/excl.mat','file')
    error("Exclusion list missing from /out/ directory! See Group_A.");
end
if ~exist('../out/pars.mat','file')
    error("File containing parameters from first study missing from /out/ directory! See Group_B.");
end
if ~exist('../out/sc_dim.mat','file')
    error("Psychiatric dimension scores missing from /out/ directory! Run E_0_process_questionnaires.m !");
end
if ~exist('../out/repeatability.mat','file')
    error('Repeatability factors not found! Run Group_C code.')
end

clear all

load ../out/pars.mat            % loads pars
load ../out/excl.mat            % loads excl
load ../out/excl_by_cond.mat    % loads excl_by_cond
load ../out/sc_dim.mat          % loads sc_dim
load ../out/repeatability.mat   % loads repeatability

parstr   = {'alpha','delta','zeta','tau'};
condstr  = {'REF','VOL','UNP'};
parrgb   = [80 80 80; 94 69 43; 46 78 93; 81 50 81]/100;
colorrgb = [41 52 29; 19 46 55; 58 22 22]/100;
facrgb   = [230 223 204; 234 191 159; 182 137 115]/255;

isout  = 1:size(pars,1);
i_excl = [];
i_excl = [i_excl excl];
for icond = 1:3
    i_excl = [i_excl excl_by_cond{icond}'];
end
isout = ismember(isout,i_excl) | isnan(sc_dim(:,1))';

% run PCA
allpars     = [pars(:,:,1) pars(:,:,2) pars(:,:,3)];
isout_pca   = ismember(isout,i_excl) | isnan(allpars(:,1))';
[cpca,spca,lpca,~,epca] = pca(zscore(allpars(~isout_pca,:),[],1));

% reform matrices to absolute subject numbering for the individual scores
for i = find(isout_pca == 1)
    if i == 1
        spca = [nan(1,size(spca,2)); spca];
    else
        spca = [spca(1:i-1,:); nan(1,size(spca,2)); spca(i:end,:)];
    end
end

npc     = 2;
rhos    = nan(3,npc);
ps      = nan(3,npc);
pns     = nan(3,npc);
ss      = cell(3,2);

% 1/ correlate psychiatric dimension scores to principal components
facstr  = {'Anxious-depressive','Compulsive behavior and intrusive thought','Social withdrawal'};
figure
clf
sgtitle(sprintf('Correlation (Psychiatric dimensions scores to PCs)\nN = %d',sum(~isout_pca)));
for idim = 1:2
    x = sc_dim(~isout,idim);
    xrange = -2:.01:4;
    for ipc = 1:npc
        subplot(2,npc,npc*(idim-1)+ipc);
        hold on
        y = spca(~isout,ipc);
        % correlation
        [pn,s] = polyfit(x,y,1);
        [py,d] = polyconf(pn,xrange,s,'alpha',0.05);
        [rho,p]  = corr(x,y,'type','Pearson');
        rhos(idim,ipc)  = rho;
        ps(idim,ipc)    = p;
        pns(idim,ipc)   = pn(1);
        ss{idim,ipc} = sigstars(p);

        shadedErrorBar(xrange,py,d,'patchSaturation',.2,'lineprops',{'LineWidth',2,'Color',facrgb(idim,:)});
        scatter(x,y,25,facrgb(idim,:),'filled');
        xline(0,':');
        yline(0,':');
        
        % calculate corrected rho for the correlation that is statistically significant
        if ipc == 1 && idim == 2
            load ../out/repeatability.mat
            r_1 = repeatability.pc.rho(1);
            n_1 = repeatability.n;
            rho_cor = rho*sqrt(1+((1-r_1)/(r_1*n_1))); % assuming no measurement noise on psychiatric dimension scores
            title(sprintf('Dimension: %s\nPC%d\nrho: %.03f; corrected rho:%.03f\np=%.03f\nslope=%.02f',facstr{idim},ipc,rho,rho_cor,p,pn(1)),'FontSize',12)
            
            % calculate corrected pearson rho using various estimations of psychiatric
            %   dimension measurement noise
            n_2 = 1000;
            fprintf('Given n_2 = %d, (number of people who retake the questionnaires)\n',n_2)
            fprintf('and r_2 is the repeatability/reliability/intraclass correlation of the psychiatric dimension,\n')
            for r_2 = .2:.1:1
                fprintf('if r_2 = %.02f, ',r_2)
                a_psy = 1+((1-r_2)/(r_2*n_2));
                rho_cor = rho*sqrt((1+((1-r_1)/(r_1*n_1))) * a_psy);
                fprintf('rho_corr = %.03f\n',rho_cor);
            end
            fprintf('Original rho = %.03f\n',rho);
        else
            title(sprintf('Dimension: %s\nPC%d\nrho: %.02f;  p=%.03f\nslope=%.02f',facstr{idim},ipc,rho,p,pn(1)),'FontSize',12)
        end
    end
    if ipc == 1
        ylabel('PC score');
    end
end

% 2/ stats on rho of correlation from psychiatric dimension score to PCs 
nq      = 2;
nb      = 1000;
pn_bs   = nan(nb,nq,npc);
rho_bs  = nan(nb,nq,npc);

% bootstrap for statistics
fprintf('Bootstrapping correlation...\n')
for ib = 1:nb
    idx = randsample(find(~isout_pca == 1),sum(~isout_pca),true);
    for idim = 1:2
        for ipc = 1:npc
            x_bs = sc_dim(idx,idim);
            y_bs = spca(idx,ipc);
            [pn_b,s] = polyfit(x_bs,y_bs,1);
            [py,d] = polyconf(pn,xrange,s,'alpha',0.05);
            [rho_bs(ib,idim,ipc),p_bs]  = corr(x_bs,y_bs,'type','Spearman');
            pn_bs(ib,idim,ipc) = pn_b(1);
        end
    end
end

p_bs    = nan(nq,2); % psy. dim, pcs
ss_ft   = cell(nq,2);
pairs   = nchoosek(1:nq,2);
nsubj   = numel(x);

% Fisher's z-transform method
for ipc = 1:npc
    ctr = 0;
    for pair = pairs'
        ctr = ctr+1;
        fprintf('Testing correlation differences between PCs on psychiatric symptoms...\n')
        z1 = atan(rhos(pair(1),ipc));
        z2 = atan(rhos(pair(2),ipc));
        z  = (z1-z2)/sqrt((nsubj-3)^-1+(nsubj-3)^-1);
        p  = normpdf(abs(z),0,1);
        ss_ft{ctr,ipc} = sigstars(p); 
        if strcmpi(ss_ft{ctr,ipc},'n.s.')
            ss_ft{ctr,ipc} = '';
        end
        p_bs(ctr,ipc) = p;
    end
end

% plot
facstr = {'AD','CIT'};
figure
clf
sgtitle(sprintf('Correlation coefficients\n(Psychiatric dimension on PCi)'))
for ipc = 1:npc
    subplot(1,2,ipc);
    hold on
    b = bar(rhos(:,ipc),'EdgeColor','none');
    b.FaceColor = 'flat';
    yline(0);
    
    for pair = pairs'
        for ctr = 1:2
            b.CData(ctr,:) = facrgb(ctr,:);
            
            % error bars from bootstrapping
            errorbar(ctr,mean(rho_bs(:,ctr,ipc),1),std(rho_bs(:,ctr,ipc),0,1),'CapSize',0,'LineWidth',2,'Color',facrgb(ctr,:)*.9);
            y_ext = (abs(mean(rho_bs(:,ctr,ipc),1))+std(rho_bs(:,ctr,ipc),0,1))*sign(rhos(ctr,ipc));
            % significance
            if rhos(ctr,ipc) < 0
                va = 'top';
            else
                va = 'bottom';
            end
            
            if strcmpi(ss(ctr,ipc),'n.s.')
                text(ctr,y_ext+.02*sign(y_ext),ss(ctr,ipc),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment',va);
            else
                text(ctr,y_ext*1.01,ss(ctr,ipc),'FontSize',20,'HorizontalAlignment','center');
            end
            
            if numel(ss_ft{ctr,ipc}) > 0 || ~isempty(ss_ft{ctr,ipc})
                plot(pair'+[.1 -.1],ones(1,2)*.38,'k','LineWidth',1.5);
            end
            
            text(mean(pair'+[.1 -.1]),.38,ss_ft(ctr,ipc),'FontSize',12,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end
    
    ylim([-.2 .42])
    if ipc == 1
        ylabel('rho')
    end
    title(sprintf('PC%d',ipc))
    set(gca,'FontSize',12)
    set(gca,'TickDir','out');
    xticks(1:2)
    xticklabels(facstr)
    xlim([.5 2.5]);
end

% 3/ stats on slopes of correlation from psychiatric dimension score to PCs 
figure
clf
sgtitle(sprintf('Slopes of correlations\n(Psychiatric dimension on PCi)'));
for idim = 1:2
    subplot(1,2,idim)
    hold on
    title(sprintf('Dimension: %s',facstr{idim}))
    b = bar(pns(idim,:),'EdgeColor','none');
    b.FaceColor      = 'flat';
    yline(0);
    y_ext = 0;
    for ipc = 1:2
        b.CData(ipc,:) = facrgb(idim,:);
        e = errorbar(ipc,mean(pn_bs(:,idim,ipc),1),std(pn_bs(:,idim,ipc),0,1),'CapSize',0,'LineWidth',2,'Color',facrgb(idim,:)*.9);
        y_ext = max(y_ext,b.YData(ipc)+e.YPositiveDelta);
        
        p = min(mean(pn_bs(:,idim,ipc)<=0),mean(pn_bs(:,idim,ipc)>=0));
        ss = sigstars(p);
        if strcmpi(ss,'n.s.')
            text(ipc,b.YData(ipc)+e.YPositiveDelta+.02,ss,'FontSize',12,'HorizontalAlignment','center');
        else
            text(ipc,b.YData(ipc)+e.YPositiveDelta+.02,ss,'FontSize',20,'HorizontalAlignment','center');
        end
    end
    p = mean((pn_bs(:,idim,1) - pn_bs(:,idim,2)) < 0);
    
    xticks(1:2)
    xlabel('principal component i');
    if idim == 1
        ylabel('slope');
    end
    set(gca,'FontSize',12);
    set(gca,'TickDir','out');
    if idim == 1
        ylim([-.35 .2])
    end 
end

% local functions
function ss = sigstars(p)
    if p <= .001
        ss = '***';
    elseif p <= .01
        ss = '**';
    elseif p <= .05
        ss = '*';
    else
        ss = 'n.s.';
    end
end
 