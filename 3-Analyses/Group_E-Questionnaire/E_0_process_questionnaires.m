% process_questionnaires

% Type:         Script
% Level:        0
% Group:        E
% Usage:        Converts scores from the mental health questionnaires into
%                   psychiatric dimension scores.
% Study:        RLVOLUNP
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

clear all

if ~exist('ques_struct.mat','file')
    error("Structure containing questionnaire scores is missing from working directory!");
end
if ~exist('factor_analysis_eigenvalues.csv','file')
    error("Factor analysis eigenvalues file missing from working directory!");
end
if ~exist('ME_phase2_excludqnadata_all.mat','file')
    error("Previous questionnaire score structure is missing from working directory!");
end

load('ques_struct.mat'); % loads ques_struct

% factor_analysis_eigenvalues.csv extracted from
%       Rouault, M., Seow, T., Gillan, C. M., & Fleming, S. M. (2018). 
%       Psychiatric symptom dimensions are associated with dissociable 
%       shifts in metacognition but not task performance. 
%       Biological psychiatry, 84(6), 443-451.
fprintf('Organizing factor analysis eigenvalues...\n');
filename = ls(fullfile('./factor_analysis_eigenvalues.csv')); 
% ------------------------------------------------------------------------
fulltable   = readtable(filename(1:end-1));
eigens      = table2array(fulltable(:,2:end));          % array of eigenvalues
[quesstr,~] = strtok(table2array(fulltable(:,1)),'.');  % array of questionnaire names

% convert questionnaire names to fit my table
%   zung -> depress
%   leb  -> social
for qstr = {'zung','leb'}
    idx = strfind(quesstr,qstr);
    idx = find(not(cellfun('isempty',idx)));
    if strcmpi(qstr,'zung')
        quesstr(idx) = {'depress'};
    elseif strcmpi(qstr,'leb')
        quesstr(idx) = {'social'};
    end
end

% extract the z-score distribution of the log transformed scores in the study of Rouault et al. (2018)
load('ME_phase2_excludqnadata_all.mat'); % loads allqna
fprintf('Calculating z-score distribution of log-transformed previous scores...\n');

labelstrorig = {'zung','anxiety','ocir','leb','bis','schizo','alcohol','eat','apathy'};
allsc_orig   = struct;

% extract data
for i = 1:numel(allqna)
    for lstr = labelstrorig
        field = getfield(allqna{i},lstr{:});
        if strcmpi(lstr{:},'leb')
            sc = field.raw.avg';
        else
            sc = field.raw';
        end
        if contains(lstr{:},{'ocir','leb','schizo','eat','alcohol'})
            sc = sc+1;
        end
        sc = log(sc);
        if i == 1
            allsc_orig = setfield(allsc_orig,lstr{:},nan(numel(sc),numel(allqna)));
        end
        allsc_orig.(lstr{:})(:,i) = sc;
    end
end
zparams	= struct;
for lstr = labelstrorig
    zparams = setfield(zparams,lstr{:},struct);
    [~,zparams.(lstr{:}).mu,zparams.(lstr{:}).sigma] = zscore(allsc_orig.(lstr{:})(:),[],'all');
end

% convert raw score on RLVOLUNP dataset to match that of the larger dataset
fprintf('Mapping raw questionnaire scores to the dataset of Rouault et al. (2018)...\n');
labelstr = {'depress','anxiety','ocir','social','bis','schizo','alcohol','eat','apathy'}; % the order of scores

% ocir,leb,shizo,eat,alcohol have +1 in the scores
% these are log transformed, and then zscored in the original dataset

nsubj   = 200;
allsc   = nan(209,nsubj);
allsc_z = nan(209,nsubj);
skip  = false;
for isubj = 1:nsubj
    if isempty(ques_struct{isubj})
        continue
    end
    for lstr = labelstr
        if ~isfield(ques_struct{isubj},lstr{:})
            skip = true;
        end
    end
    if skip
        skip = false;
        continue
    end
    for lstr = labelstr
        ind = zeros(209,1);
        idx = strfind(quesstr,lstr);
        idx = find(not(cellfun('isempty',idx)));
        ind(idx) = 1;
        field = getfield(ques_struct{isubj},lstr{:});
        if strcmpi(lstr{:},'social')
            sc = field.avg';
        else
            sc = field.raw';
        end
        if contains(lstr{:},{'ocir','social','schizo','eat','alcohol'})
            sc = sc+1;
        end 
        sc = log(sc);
        allsc(logical(ind),isubj) = sc;
        if strcmpi(lstr{:},'depress')
            zmu  = zparams.('zung').mu;
            zsig = zparams.('zung').sigma;
        elseif strcmpi(lstr{:},'social')
            zmu  = zparams.('leb').mu;
            zsig = zparams.('leb').sigma;
        else
            zmu  = zparams.(lstr{:}).mu;
            zsig = zparams.(lstr{:}).sigma;
        end
        % z-score based on large dataset
        allsc_z(logical(ind),isubj) = (sc - zmu)/zsig;
    end
end

fprintf('Calculating scores on psychiatric dimension...\n');
% score on psychiatric dimension 
sc_dim = nan(nsubj,3);
for idim = 1:3
    sc_dim(:,idim) = sum(eigens(:,idim).*allsc_z,1);
end

% save psychiatric dimension scores
save('../out/sc_dim.mat','sc_dim');
