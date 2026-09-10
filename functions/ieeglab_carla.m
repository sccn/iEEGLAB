function [Vout, CAR, stats] = ieeglab_carla(tt, V, srate, opts)
% ieeglab_carla() - Common Average Referencing by Least Anticorrelation (CARLA).
%
% Channels are ranked by increasing mean cross-trial covariance (variance if a
% single trial). Channels are then iteratively added to a candidate common
% average, and at each subset size the mean correlation between each
% un-re-referenced channel and all re-referenced channels is computed on
% bootstrapped mean signals. zMin is the mean correlation belonging, at each
% bootstrap, to the most anticorrelated ("most responsive") channel. The optimal
% CAR size is the subset size at which anticorrelation is least, i.e. where
% zMinMean is at its maximum (least negative).
%
% The intuition: too few channels in the common average over-represents the
% noise of each constituent; too many pulls genuinely responsive channels into
% the reference. The peak in the test statistic finds the largest subset that
% stops short of the responsive channels.
%
% Usage:
%   [Vout, CAR, stats] = ieeglab_carla(tt, V, srate)
%   [Vout, CAR, stats] = ieeglab_carla(tt, V, srate, opts)
%
% Inputs:
%   tt    - [1 x t] time points matching V, in SECONDS
%   V     - [n x t x k] signal, n channels by t time points by k trials
%           ([n x t] accepted when k == 1). Channels that must not enter the
%           common average (stimulated contacts, bad channels) should be set to
%           NaN across all trials; they are excluded from ranking and from the
%           reference, but are still re-referenced on output.
%   srate - sampling frequency in Hz. Strongly recommended; if omitted it is
%           estimated from tt.
%   opts  - (optional) struct:
%      .winResp    [1x2] response window in s for ranking/correlation. Default [0.01 0.3]
%      .nboot      integer, bootstrap samples. Default 100. Ignored if k == 1.
%      .sens       logical, use the more sensitive cutoff (before the first
%                  significant decrease in zMinMean) instead of the global
%                  maximum. Default false.
%      .notchFirst logical, notch the signal before ranking only (the returned
%                  data is re-referenced from the UNfiltered input). Default true.
%      .lineFreq   line noise fundamental in Hz, 50 or 60. Default 60.
%      .minCARsize integer, smallest allowed common average. Default 2.
%      .verbose    logical. Default true.
%
% Outputs:
%   Vout  - re-referenced signal, same shape as V
%   CAR   - [t x k] the common average subtracted at each trial
%   stats - struct with fields chsUsed, vars, order, zMin, zMinMean, nOptimum
%
% Adapted for iEEGLAB from CARLA.m (2023/04/27 by Harvey Huang), part of the
% CARLA manuscript package, https://github.com/hharveygit/CARLA_JNM
%
% If you use this function, please cite:
%   Huang H, Ojeda Valencia G, Gregg NM, Osman GM, Montoya MN, Worrell GA,
%   Miller KJ, Hermes D (2024). CARLA: Adjusted common average referencing for
%   cortico-cortical evoked potential data. Journal of Neuroscience Methods,
%   407:110153. https://doi.org/10.1016/j.jneumeth.2024.110153
%
% Copyright (C) 2023 Harvey Huang
% Adapted for iEEGLAB by Cedric Cannard, 2026
%
% This program is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version. See <https://www.gnu.org/licenses/>.

% ---------------- options ----------------
if nargin < 4, opts = struct(); end
def = struct('winResp',[0.01 0.3], 'nboot',100, 'sens',false, ...
             'notchFirst',true, 'lineFreq',60, 'minCARsize',2, 'verbose',true);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opts, fn{ii}) || isempty(opts.(fn{ii}))
        opts.(fn{ii}) = def.(fn{ii});
    end
end

% Accept [t x n] when k == 1, matching the reference implementation
if ismatrix(V) && size(V,1) == length(tt) && size(V,2) ~= length(tt)
    V = V';
end

tt = tt(:)';
nChs = size(V,1);
nTrs = size(V,3);
assert(size(V,2) == length(tt), ...
    'ieeglab_carla:badDims', ...
    'Second dimension of V (%d) does not match tt (%d). Expected channels x time x trials.', ...
    size(V,2), length(tt));

if nargin < 3 || isempty(srate)
    srate = (length(tt)-1) / (max(tt)-min(tt));
    warning('ieeglab_carla:estimatedSrate', ...
        'Sampling frequency estimated as %.2f Hz from tt. Pass srate explicitly.', srate);
end

stats = struct();

% Channels excluded from the reference: all-NaN across the whole epoch set.
% These are still re-referenced on output, they just cannot join the average.
isExcluded = all(all(isnan(V), 2), 3);
isExcluded = isExcluded(:);
poolIdx = find(~isExcluded);
if numel(poolIdx) < opts.minCARsize
    warning('ieeglab_carla:tooFewChannels', ...
        'Only %d channel(s) available for the common average (need >= %d). Returning data unchanged.', ...
        numel(poolIdx), opts.minCARsize);
    Vout = V;
    CAR  = zeros(size(V,2), nTrs);
    stats.chsUsed = []; stats.vars = nan(nChs,1); stats.order = poolIdx;
    stats.zMinMean = []; stats.zMin = []; stats.nOptimum = 0;
    return
end

% ---------------- notch filter for ranking only ----------------
Vclean = V;
if opts.notchFirst
    nyq = srate/2;
    haveDSP = exist('designfilt','file') == 2 && license('test','Signal_Toolbox');
    if ~haveDSP
        if opts.verbose
            warning('ieeglab_carla:noDSP', ...
                'Signal Processing Toolbox not available; skipping the notch used for channel ranking.');
        end
    else
        for ff = opts.lineFreq : opts.lineFreq : 3*opts.lineFreq
            if ff+2 >= nyq, break; end   % harmonic above Nyquist
            try
                dNotch = designfilt('bandstopiir', 'FilterOrder', 4, ...
                    'DesignMethod', 'butter', ...
                    'HalfPowerFrequency1', ff-2, ...
                    'HalfPowerFrequency2', ff+2, ...
                    'SampleRate', srate);
                for ii = 1:nTrs
                    trial = Vclean(:,:,ii)';
                    ok = ~any(isnan(trial), 1);     % filtfilt propagates NaN across the trial
                    trial(:,ok) = filtfilt(dNotch, trial(:,ok));
                    Vclean(:,:,ii) = trial';
                end
            catch ME
                if opts.verbose
                    warning('ieeglab_carla:notchFailed', ...
                        'Notch at %g Hz failed (%s); continuing without it.', ff, ME.message);
                end
                break
            end
        end
    end
end

% ---------------- response segment ----------------
tmask = tt >= opts.winResp(1) & tt <= opts.winResp(2);
if nnz(tmask) < 3
    error('ieeglab_carla:emptyWindow', ...
        'Response window [%g %g] s selects %d sample(s) of tt (range [%g %g] s). Widen it or check units.', ...
        opts.winResp(1), opts.winResp(2), nnz(tmask), tt(1), tt(end));
end
Vseg = Vclean(:, tmask, :);

% ---------------- rank channels ----------------
% Single trial -> variance. Multiple trials -> mean cross-trial covariance,
% which is a better estimator of consistent (i.e. evoked) structure.
stats.vars = nan(nChs,1);
if nTrs == 1
    stats.vars(poolIdx) = var(Vseg(poolIdx,:), 0, 2);
else
    for ii = poolIdx'
        covCurr = cov(squeeze(Vseg(ii,:,:)));
        stats.vars(ii) = mean(covCurr(logical(triu(ones(size(covCurr)),1))), 'all');
    end
end
[~, ordLocal] = sort(stats.vars(poolIdx));
stats.order = poolIdx(ordLocal);         % global channel indices, increasing covariance
nPool = numel(stats.order);

% ---------------- sweep CAR size ----------------
if nTrs == 1
    stats.zMin = nan(nPool, nPool);
else
    stats.zMin = nan(nPool, nPool, opts.nboot);
end

for ii = max(2, opts.minCARsize) : nPool

    sub = stats.order(1:ii);
    VsegReref = Vseg - mean(Vseg(sub,:,:), 1);

    if nTrs == 1
        Useg      = Vseg(sub,:)';
        UsegReref = VsegReref(sub,:)';
        r = corr(Useg, UsegReref);
        r(1:(ii+1):end) = nan;               % drop self-self
        z = atanh(r);
        [~, kkMost] = min(mean(z, 2, 'omitnan'));
        stats.zMin(1:ii, ii) = z(kkMost,:)';
        continue
    end

    for jj = 1:opts.nboot
        inds = randi(nTrs, 1, nTrs);          % bootstrap resample of trials
        Useg      = mean(Vseg(sub,:,inds), 3)';
        UsegReref = mean(VsegReref(sub,:,inds), 3)';
        r = corr(Useg, UsegReref);
        r(1:(ii+1):end) = nan;
        z = atanh(r);
        [~, kkMost] = min(mean(z, 2, 'omitnan'));
        stats.zMin(1:ii, ii, jj) = z(kkMost,:)';
    end
end

stats.zMinMean = mean(stats.zMin, 1, 'omitnan');   % 1 x nPool x nboot

% ---------------- pick the optimum ----------------
if opts.sens && nTrs > 1
    nOptimum = local_sensitive_optimum(stats, nPool, opts);
else
    [~, nOptimum] = max(mean(stats.zMinMean, 3));
end
nOptimum = max(nOptimum, min(opts.minCARsize, nPool));

stats.zMinMean = squeeze(stats.zMinMean);
stats.nOptimum = nOptimum;
stats.chsUsed  = sort(stats.order(1:nOptimum));

% ---------------- apply, on the UNfiltered input ----------------
CAR  = mean(V(stats.chsUsed,:,:), 1, 'omitnan');
Vout = V - CAR;
CAR  = squeeze(CAR);

if opts.verbose
    fprintf('[CARLA] %d/%d channels in the common average (ranked by %s).\n', ...
        nOptimum, nPool, ternary(nTrs==1,'variance','cross-trial covariance'));
end

end

% ================= local helpers =================

function nOptimum = local_sensitive_optimum(stats, nPool, opts)
% Optimum at the last peak before the first STATISTICALLY significant decrease
% in zMinMean, rather than the global maximum. More sensitive to small numbers
% of responsive channels.
nMin = max(ceil(0.1*nPool), 2);
zMMxTrs = mean(stats.zMinMean(1,:,:), 3);
ii = nMin;
nOptimum = nPool;
while ii <= nPool
    if ii == nPool, nOptimum = ii; break; end
    if zMMxTrs(ii+1) > zMMxTrs(ii), ii = ii + 1; continue; end

    zMMxTrs(1:ii-1) = nan;
    nextGreater = find(zMMxTrs > zMMxTrs(ii), 1, 'first');
    if isempty(nextGreater), nOptimum = ii; break; end

    zMMCurr   = squeeze(stats.zMinMean(1, ii, :));
    [~, idx]  = min(zMMxTrs(1:nextGreater-1));
    zMMTrough = squeeze(stats.zMinMean(1, idx, :));

    % Left-tailed 95% CI on the pairwise difference (trough - peak). The
    % bootstrap samples are not independent, so we test the difference of
    % zMinMean directly rather than the mean of zMinMean.
    [X, Y] = meshgrid(zMMTrough, zMMCurr);
    diffs = X - Y;
    if prctile(diffs(:), 95) < 0
        nOptimum = ii; break
    end
    ii = nextGreater;
end
if nOptimum == nMin
    warning('ieeglab_carla:floorReached', ...
        'Optimum CAR detected at the 10%% minimum-size floor (%d channels).', nMin);
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
