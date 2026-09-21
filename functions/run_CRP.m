function [crp_parms, crp_projs] = run_CRP(V, t_win, opts)
% RUN_CRP  Canonical Response Parameterization (robust, self-contained)
%
%   [crp_parms, crp_projs] = RUN_CRP(V, t_win, opts)
%
%   Perform Canonical Response Parameterization (CRP) on single-pulse,
%   stimulation-evoked voltage responses across trials. The method computes
%   cross-projection profiles over increasing response durations, selects an
%   optimal duration, and parameterizes trials using the first principal
%   component ("canonical shape"). Significance metrics are provided from the
%   projection distributions.
%
%   INPUTS
%     V      : Single-trial voltage matrix.
%              Accepted orientations:
%                • [T x K]  time points × trials  (preferred)
%                • [K x T]  trials × time points  (auto-transposed)
%              (Rows must correspond to the time vector passed in t_win.)
%
%     t_win  : Time vector (seconds) of length T (row or column OK),
%              corresponding to the time axis of V.
% 
%   opts   : (optional) struct with fields:
%              .verbose (default=true)  - print progress
%              .label   (default='')    - label to include in headers
% 
%   OUTPUTS
%     crp_projs : struct with projection profiles & stats
%        .proj_tpts         — 1×M time points (seconds) at which projections
%                              were computed (grows from 10 samples by t_step)
%        .S_all             — (K^2−K)×M set of cross-projection magnitudes
%        .mean_proj_profile — 1×M mean projection profile over durations
%        .var_proj_profile  — 1×M variance of projection profile
%        .tR_index          — scalar index into proj_tpts of chosen duration
%        .avg_trace_input   — T×1 simple average trace over the full window
%        .stat_indices      — indices used for non-overlapping statistics
%        .t_value_tR        — t-statistic at the chosen duration τ_R
%        .p_value_tR        — one-tailed p-value at τ_R
%        .t_value_full      — t-statistic at the full window
%        .p_value_full      — one-tailed p-value at the full window
%
%     crp_parms : struct with canonical parameterization
%        .V_tR         — tR×K reduced matrix (truncated to τ_R)
%        .al           — 1×K alpha weights projecting V_tR onto C
%        .C            — tR×1 canonical shape (1st PC of V_tR)
%        .ep           — tR×K residual after removing C*al
%        .tR           — scalar, response duration (seconds)
%        .parms_times  — 1×tR time vector for parameterized data
%        .avg_trace_tR — tR×1 average trace (truncated to τ_R)
%        .al_p         — 1×K alpha-prime (alpha / sqrt(#samples))
%        .epep_root    — 1×K residual norm per trial
%        .Vsnr         — 1×K signal-to-noise per trial
%        .expl_var     — 1×K explained variance per trial
%        .cod          — 1×K coefficient of determination per trial
%
%   ASSUMPTIONS & REQUIREMENTS
%   • V should be baseline-corrected and filtered upstream.
%   • Minimal data: ≥10 time samples and ≥2 trials; otherwise function returns
%     well-formed, empty outputs (no crash).
%   • Sampling rate is inferred from t_win (assumed ~uniform spacing).
%   • The projection step scales to sqrt(seconds) for unit consistency.
%
%   USAGE (pipeline context)
%       % sig: [chan x time x trials], timevec: 1×T
%       pick   = (timevec >= 0.015) & (timevec <= 1.0);  % CRP window
%       V_full = squeeze(sig(ii, :, :));   % [time x trials] for channel ii
%       V_win  = V_full(pick, :);          % time-windowed
%       [crp_parms, crp_projs] = run_CRP(V_win, timevec(pick));
%
%   IMPLEMENTATION NOTES
%   • Orientation-agnostic: the function transposes V if needed so that
%     size(V,1) == numel(t_win) before processing.
%   • Projection durations are sampled at:
%       proj_tpts = 10:t_step:T        (default t_step = 5 samples)
%     You can tune t_step inside this file for smoother profiles vs. speed.
%   • Linear "kernel-trick" PCA (kt_pca) is used to obtain the canonical shape.
%   • If Statistics Toolbox is missing, p-values are computed via tcdf fallback.
%
%     Original method & implementation: Kai J. Miller, M.D., Ph.D. (2022; revised 2023)
%
%     Adaptation for iEEGLAB plugin: Cedric Cannard (2025)
%
%   REFERENCE: 
%       Miller, K. J., et al. (2023). Canonical Response Parameterization: 
%       Quantifying the structure of responses to single-pulse intracranial 
%       electrical brain stimulation. PLoS computational biology, 19(5), 
%       e1011105.
% 
%   LICENSE
%   Copyright (C) 2022 Kai J. Miller
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License, version 3 or later.
%   See <https://www.gnu.org/licenses/> for details.


% ---------- options & verbosity ----------
if nargin < 3, opts = struct(); end
opts = local_merge_opts(opts, struct( ...
    'verbose',    true, ...
    't_step',     5, ...
    'min_T',      10, ...
    'min_K',      2, ...
    'nan_policy', 'omittrial' ...
));

vprint = @(varargin) (opts.verbose && fprintf(varargin{:}));

% Header respects the verbose flag: this function is called once per
% channel-by-site pair, so an unconditional banner floods the command window
% during a scripted run.
if opts.verbose
    if isfield(opts,'label')
        fprintf('\n================= [CRP] Run: %s =================\n', opts.label);
    else
        fprintf('\n================= [CRP] New run =================\n');
    end
end

% ---------- Input hygiene ----------
if isrow(t_win), t_win = t_win(:); end
T1 = size(V,1); T2 = size(V,2);
if T1 == numel(t_win)
    % V already [T x K]
elseif T2 == numel(t_win)
    V = V.';  % transpose to [T x K]
else
    error('run_CRP:badDims','Neither dimension of V matches numel(t_win).');
end
[T,K] = size(V);

vprint('[CRP] Oriented V to [time x trials]: T=%d, K=%d\n', T, K);

% ---------- Basic checks ----------
if any(diff(t_win) <= 0)
    error('run_CRP:badTime','t_win must be strictly increasing.');
end
dt = diff(t_win);
srate = 1/mean(dt);
rel_jitter = std(dt)/mean(dt);
if rel_jitter > 1e-3
    vprint('[CRP] WARNING: non-uniform sampling (std(dt)/mean(dt)=%.2e). Proceeding.\n', rel_jitter);
end

% NaN handling per policy
if any(~isfinite(V(:)))
    switch lower(opts.nan_policy)
        case 'omittrial'
            keep = all(isfinite(V),1);
            dropped = find(~keep);
            if ~isempty(dropped)
                vprint('[CRP] Omitting %d/%d trial(s) due to NaNs/Infs.\n', numel(dropped), K);
            end
            V = V(:, keep);
        case 'zero'
            vprint('[CRP] Replacing %d NaN/Inf entries with 0.\n', sum(~isfinite(V(:))));
            V(~isfinite(V)) = 0;
        case 'error'
            error('run_CRP:NaN','V contains NaN/Inf and nan_policy="error".');
        otherwise
            error('run_CRP:badOption','Unknown nan_policy: %s', opts.nan_policy);
    end
    [T,K] = size(V); % update after omission
end

% Minimal data checks
if T < opts.min_T || K < opts.min_K || any(isnan(t_win))
    vprint('[CRP] Insufficient data (T=%d, K=%d). Returning empty outputs.\n', T, K);
    [crp_parms, crp_projs] = crp_empty_outputs(t_win, T, K);
    return
end

% Small DC/baseline drift warning (assumes upstream baseline removal)
mcol = mean(V,1);
if any(abs(mcol) > 1e-6) % unit-agnostic tiny threshold
    vprint('[CRP] Note: non-zero per-trial mean detected (baseline already removed upstream?).\n');
end

% ---------- Projection profiles ----------
t_step    = max(1, round(opts.t_step));    % force sane integer
proj_tpts = 10:t_step:T;                   % sample indices for profiles
if isempty(proj_tpts), proj_tpts = 10:T; end
M = numel(proj_tpts);

vprint('[CRP] Computing projection profiles: t_step=%d samples, M=%d durations (srate≈%.2f Hz)\n', t_step, M, srate);

m    = zeros(1,M);
v2   = zeros(1,M);
S_all = NaN(K^2 - K, M);

% lightweight progress updates
tick = max(1, floor(M/10));
for j = 1:M
    k = proj_tpts(j);
    S = ccep_proj(V(1:k,:));       % 1 x (K^2-K)
    S = S / sqrt(srate);           % -> sqrt(seconds)
    S_all(:,j) = S(:);
    m(j)  = mean(S);
    v2(j) = var(S);
    % if opts.verbose && (mod(j,tick)==0 || j==M)
    %     fprintf('[CRP]  ...duration %3d/%3d (k=%d samples)\n', j, M, k);
    % end
end
[~,tt] = max(m);
tR_samp = proj_tpts(tt);
tR_sec  = t_win(tR_samp);

vprint('[CRP] Selected response duration: tR=%d samples (%.3f s)\n', tR_samp, tR_sec);

% ---------- Parameterize at response duration ----------
V_tR     = V(1:tR_samp,:);           % [tR x K]
[E_tR,~] = kt_pca(V_tR);             % linear kernel PCA
C        = E_tR(:,1);                % canonical shape (unit-normalized)
al       = C.' * V_tR;               % 1 x K
ep       = V_tR - C * al;            % tR x K

% ---------- Outputs ----------
% Projections / stats
crp_projs.proj_tpts         = t_win(proj_tpts);
crp_projs.S_all             = S_all;
crp_projs.mean_proj_profile = m;
crp_projs.var_proj_profile  = v2;
crp_projs.tR_index          = tt;
crp_projs.avg_trace_input   = mean(V,2);
crp_projs.stat_indices      = get_stat_indices(K);
crp_projs.t_value_tR        = tstat_safe(S_all(crp_projs.stat_indices,tt));
crp_projs.p_value_tR        = pval_safe(S_all(crp_projs.stat_indices,tt));
crp_projs.t_value_full      = tstat_safe(S_all(crp_projs.stat_indices,end));
crp_projs.p_value_full      = pval_safe(S_all(crp_projs.stat_indices,end));
% Useful metadata
crp_projs.settings   = struct('t_step',t_step,'srate',srate,'min_T',opts.min_T,'min_K',opts.min_K,'nan_policy',opts.nan_policy);

% Parameterizations
crp_parms.V_tR         = V_tR;
crp_parms.al           = al;
crp_parms.C            = C;
crp_parms.ep           = ep;
crp_parms.tR           = tR_sec;
crp_parms.parms_times  = t_win(1:tR_samp);
crp_parms.avg_trace_tR = mean(V_tR,2);
crp_parms.al_p         = al / sqrt(length(C));
crp_parms.epep_root    = sqrt(diag(ep.'*ep)).';
crp_parms.Vsnr         = al ./ sqrt(diag(ep.'*ep)).';
crp_parms.expl_var     = 1 - (diag(ep.'*ep)).' ./ (diag(V_tR.'*V_tR)).';
crp_parms.cod          = 1 - (diag(ep.'*ep)).' ./ (diag( (V_tR-mean(V_tR)).'*(V_tR-mean(V_tR)) ).');

% Quick summary
vprint('[CRP] Summary: K=%d trials, mean explained var=%.3f, tR=%.3f s, t(full)=%.3f (p=%.3g)\n', ...
    K, mean(crp_parms.expl_var,'omitnan'), crp_parms.tR, crp_projs.t_value_full, crp_projs.p_value_full);

if opts.verbose
    fprintf('-------------------------------------------------------------------\n\n');
end

% ================== local helpers ==================
function S0 = ccep_proj(Vin)
    % normalize each trial (L2), handle divide-by-zero
    Vin0 = Vin ./ (ones(size(Vin,1),1) * (sum(Vin.^2,1).^.5));
    Vin0(~isfinite(Vin0)) = 0;
    P = Vin0.' * Vin;                   % semi-normalized projections
    P(1:size(P,1)+1:end) = NaN;         % remove diagonal
    S0 = reshape(P,1,[]); S0(isnan(S0)) = [];

function [E,S] = kt_pca(X)
    % X: [T x K] (time x trials)
    [F,S2] = eig(X.'*X);                                 % KxK
    [svals, idx] = sort(real(diag(S2)),'descend');
    F = F(:, idx);                                       % sorted
    S = sqrt(max(svals, 0));                             % Kx1 ≥0
    ES = X * F;                                          % [T x K]
    den = repmat(max(S, eps).', size(X,1), 1);           % [T x K]
    E  = ES ./ den;                                      % unit-normalized eigenvectors

function stat_indices = get_stat_indices(N)
    % Non-overlapping indices over off-diagonal entries (Miller 2023)
    stat_indices = 1:2:(N^2-N);
    if rem(N,2)==1
        b = zeros(size(stat_indices));
        for k = 1:N
            if mod(k,2)==0
                b(((k-1)*((N-1)/2)+1):(k*((N-1)/2))) = 1;
            end
        end
        stat_indices = stat_indices + b;
    end

function t = tstat_safe(x)
    x = x(:);
    sx = std(x,'omitnan');
    n  = sum(isfinite(x));
    if n < 2 || sx == 0, t = NaN; else, t = mean(x,'omitnan')/(sx/sqrt(n)); end

function p = pval_safe(x)
    x = x(:);
    x = x(isfinite(x));
    if numel(x) < 2, p = NaN; return; end
    if exist('ttest','file') == 2
        try, [~,p] = ttest(x,0,'tail','right'); catch, p = NaN; end
    else
        t = tstat_safe(x);
        df = max(numel(x)-1,1);
        if isfinite(t), p = 1 - tcdf(t, df); else, p = NaN; end
    end

function [parms, projs] = crp_empty_outputs(tw, Tloc, Kloc)
    if isempty(tw), tw = nan(Tloc,1); end
    projs.proj_tpts = tw(:).';
    projs.S_all = [];
    projs.mean_proj_profile = [];
    projs.var_proj_profile  = [];
    projs.tR_index = [];
    projs.avg_trace_input = [];
    projs.stat_indices = [];
    projs.t_value_tR = NaN; projs.p_value_tR = NaN;
    projs.t_value_full = NaN; projs.p_value_full = NaN;
    projs.settings = struct();
    parms.V_tR = []; parms.al = []; parms.C = []; parms.ep = [];
    parms.tR = NaN; parms.parms_times = []; parms.avg_trace_tR = [];
    parms.al_p = []; parms.epep_root = []; parms.Vsnr = [];
    parms.expl_var = []; parms.cod = [];

function opts = local_merge_opts(opts, defaults)
    f = fieldnames(defaults);
    for i = 1:numel(f)
        if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
            opts.(f{i}) = defaults.(f{i});
        end
    end