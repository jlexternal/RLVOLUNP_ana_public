% gen_drift

% Type:         Function
% Level:        0
% Group:        TaskGen
% Usage:        Generates the mean drift values for the VOL condition and the static
%                   mean values for the REF and UNP conditions
% Study:        RLVOLUNP
%
% Requirements: betafun.m, betapar.m, betasmp.m
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

function out = gen_drift(cfg)

addpath('./functions_TaskGen/');

% load configuration into local
ngen        = cfg.ngen;         % number of episodes to be generated
nepis       = cfg.nepis;        % number of episodes desired
epimin      = cfg.epimin;       % minimum episode length
epimax      = cfg.epimax;       % maximum episode length
avgmin      = cfg.avgmin;       % minimum average value allowed
avgmax      = cfg.avgmax;       % maximum average value allowed
var_drift   = cfg.var_drift;    % drifting variance

%fprintf('generating episodes\n');

get_pr = @(m,v)betasmp(m,v); % sample from beta distr. w/ mean m and variance v

pr_all      = cell(ngen,1);
pr          = zeros(1,epimax+2);

while true 
    % Generate mean-drift episodes
    for i = 1:ngen
        if mod(i,100) == 0
%            fprintf('Generating %d of %d episodes...\n',i,ngen);
        end

        while true
            pr(:) = 0;
            pr(1) = get_pr(0.5,var_drift);

            % initialize episode value
            if pr(1) < 0.5
                pr(1) = 1-pr(1);
            end
            % create random walk from initial value
            for j = 2:epimax+2
                pr(j) = get_pr(pr(j-1),var_drift);
                if pr(j) < 0.5
                    break
                end
            end

            epilen = j-1;
            epiavg = mean(pr(1:j-1));
            if epilen >= epimin && ...          % episode length is at least at minimum length
               epilen <= epimax && ...          % episode length is at most at maximum length
               epiavg >= avgmin && ...          % average value is at least at minimum value
               epiavg <= avgmax && ...          % average value is at most at maximum value
               max(pr(1:j-1)) <= 0.90 && ...    % max value within episode does not exceed preset max
               max(abs(diff(pr(1:j)))) <= 0.10  % trial to trial differences are not overly great
                break
            end
        end
        pr_all{i} = pr(1:j-1);
    end

    % Filter through episodes whose median lies within some range
    xavg = cellfun(@median,pr_all);
    i = xavg >= 0.55 & xavg <= 0.65;
    pr_all = pr_all(i);
    
    % get desired number of episodes
    if numel(pr_all) > 1000
        pr_all = pr_all(1:1000);
        break
    else
        disp('1000 episodes within criteria not found. Re-running episode generation...\n');
    end
end
pr_all_vec      = horzcat(pr_all{:});

% Identify start/switch points 
nepilen     = cellfun(@numel,pr_all);
epistart    = nan(1,nepis);
epistart(1) = 1;
for i = 2:nepis
    epistart(i) = epistart(i-1)+nepilen(i-1);
end

epis        = pr_all_vec;
epis_vol    = epis;
epis_avg    = nan(1,nepis); % average value of episode 
sw_trs      = false(1,sum(nepilen));

% invert values at switch points
for i = 1:nepis
    epis_avg(i) = mean(epis_vol(epistart(i):(epistart(i)+nepilen(i)-1)));
    
    % invert values for VOL condition
    if mod(i,2) == 0
        epis_vol(epistart(i):(epistart(i)+nepilen(i)-1)) = 1-epis_vol(epistart(i):(epistart(i)+nepilen(i)-1));
        epis_avg(i) = 1-epis_avg(i);
        
        sw_trs(epistart(i):(epistart(i)+nepilen(i)-1)) = 1;
    end
    
    % create mean "drift" for the REF condition
    epis_ref(epistart(i):(epistart(i)+nepilen(i)-1)) = epis_avg(i);
end

% Get drifting variance
var_drift_ef = var(epis_vol(2:end)-epis_vol(1:end-1));
%fprintf('True drift variance: %.4f\nEffective drift variance: %.4f\n',var_drift,var_drift_ef);

out             = struct;
out.epis_vol    = epis_vol;     % mean drift value for the VOL condition
out.epis_avg    = epis_avg;     % static values for REF and UNP conditions
out.epis_ref    = epis_ref;     % mean static values for the REF & UNP conditions
out.vd_eff      = var_drift_ef; % effective drift variance for the VOL drift values
out.epistart    = epistart;     % indices of switch/start points of the drifts
out.nepilen     = nepilen;      % lengths of episodes contained within entire drift
out.sw_trs      = sw_trs;       % switch point indices 
out.cfg_in      = cfg;

end