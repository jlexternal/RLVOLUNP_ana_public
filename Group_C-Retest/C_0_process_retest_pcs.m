% process_retest_pcs

% Type:         Script
% Level:        0
% Group:        C
% Usage:        Calculate reliability/repeatability of the principal components
%                   extracted from the initial study
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021
clear all
addpath('../out')

if ~exist('out_pca_test.mat','file')
    error("Principal component analysis results missing from /out/ directory!\nRun Group B analyses first!");
end
if ~exist('pars.mat','file')
    error("File containing parameters from first study missing from /out/ directory!");
end
if ~exist('pars_rt.mat','file')
    error("File containing parameters from retest missing from /out/ directory!");
end

load out_pca_test.mat   % loads cpca, spca
load pars.mat           % loads pars
load pars_rt.mat        % loads pars_rt
load sc_dim.mat         % loads sc_dim
load excl.mat           % loads excl
load excl_by_cond.mat   % loads excl_by_cond

i_excl = [];
i_excl = [i_excl excl];
for icond = 1:3
    i_excl = [i_excl excl_by_cond{icond}'];
end
i = setdiff(1:size(pars,1),i_excl);
spca_t = nan(size(pars,1),size(spca,2));

spca_t(i,:) = spca;
cpca_t = cpca;

clearvars spca cpca

% repeatability/reliability of PC1, PC2
parrgb      = [80 80 80; 94 69 43; 46 78 93; 81 50 81]/100;
colorrgb    = [41 52 29; 19 46 55; 58 22 22]/100;
condstr     = {'REF','VOL','UNP'};

allpars_rt = [pars_rt(:,:,1) pars_rt(:,:,2) pars_rt(:,:,3)];
idx = find(~isnan(allpars_rt(:,1)));

spca_rt = nan(200,12);
spca_rt(idx,:) = zscore(allpars_rt(idx,:),[],1)*cpca_t;

repeatability = struct;
repeatability.n = 100;
repeatability.rho = struct;

npc = 2;
for ipc = 1:npc
    figure
    x = spca_t(idx,ipc);
    y = spca_rt(idx,ipc);
    xleft   = min(x)-.1;
    xright  = max(x)+.1;
    [pn,s]  = polyfit(x,y,1);
    [py,d]  = polyconf(pn,xleft:.01:xright,s,'alpha',0.05);
    [rho,p] = corr(x,y,'Type','Spearman');
    repeatability.pc.rho(ipc) = rho;
    scatter(x,y)
    hold on
    plot(xleft:xright,xleft:xright,':','Color','k')
    shadedErrorBar(xleft:.01:xright,py,d,'patchSaturation',.1,'lineprops',{'LineWidth',2});
    xlabel('Test');
    ylabel('Retest');
    title(sprintf('PC%d scores (test-retest)\nrho: %.02f\np=%.03f\nslope=%.02f',ipc,rho,p,pn(1)));
    set(gca,'TickDir','out');
    set(gca,'FontSize',12);
end

save('../out/repeatability.mat','repeatability');