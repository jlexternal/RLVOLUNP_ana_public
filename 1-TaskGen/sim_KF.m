% sim_KF

% Type:         Function
% Level:        0
% Group:        TaskGen
% Usage:        Simulates a block of the task via a Kalman Filter
% Study:        RLVOLUNP
%
% Usage:        Match difficulty between VOL and UNP conditions (in gen_sessions.m)
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021


function out = sim_KF(cfg)

% cfg       : configuration input
% epis      : mean drift array
% traj      : sampled trajectory array
% isvol     : boolean: volatile condition or not 
% vs        : measurement uncertainty / spread around mean
% vd        : process noise / drift variance
% epistart  : indices of start/switch trial

epis        = cfg.epis;
traj        = cfg.traj;
isvol       = cfg.isvol;
epistart    = cfg.epistart;

if isfield(cfg,'vs')
    vs = cfg.vs;
else
    vs = var(epis-traj);
end
if isfield(cfg,'vd')
    vd = cfg.vd;
else
    vd = var(epis(2:end)-epis(1:end-1));
end

nt = numel(traj);

% allocate variables
kt = nan(1,nt); % kalman gain
vt = nan(1,nt); % posterior variance
mt = nan(1,nt); % estimated mean
rt = nan(1,nt); % KF responses
pt = nan(1,nt); % probability of correct KF responses (argmax)

resp_up = false; % whether response is above .5 or below

% filter
for it = 1:nt
    if ismember(it,epistart)
        resp_up = ~resp_up;
    end
    
    if it == 1 || (~isvol && epis(it) ~= epis(it-1))
        mt(it) = mean(traj); % mean of latent trajectory
        vt(it) = var(traj);  % variance of latent trajectory
        kt(it) = 1;
    else
        mt(it) = mt(it-1);
        vt(it) = vt(it-1);
        kt(it) = kt(it-1);
    end
    
    % get responses and correct response probabilities for the KF
    if resp_up
        rt(it) = mt(it) >= .5;
        pt(it) = 1-normcdf(.5,mt(it),vt(it));
    else
        rt(it) = mt(it) < .5;
        pt(it) = normcdf(.5,mt(it),vt(it));
    end
    
    % KF update
    kt(it) = vt(it)/(vt(it)+vs);
    mt(it) = mt(it) + kt(it)*(traj(it)-mt(it));
    vt(it) = (1-kt(it))*vt(it)+vd;
end

out.cfg = cfg;
out.mt  = mt;
out.vt  = vt;
out.kt  = kt;
out.rt  = rt;
out.pt  = pt;

end