function [T, R] = ieeglab_band_stats(S, varargin)
% ieeglab_band_stats() - Permutation test of band power against the baseline,
%                        for each contact (high-frequency activity by default),
%                        with correction across contacts.
%
% Usage:
%   S = ieeglab_spectral(EEG, 'type', 'timefreq', 'keeptrials', 'band');
%   [T, R] = ieeglab_band_stats(S);
%   EEG = pop_ieeglab_spectral(EEG, 'type', 'timefreq', 'stats', true);   % the same, with figures
%   [T, R] = ieeglab_band_stats(S, 'window', [100 500], 'tail', 'positive', 'per_time', true);
%
% For each contact, epoch and frequency of the band, the change is the mean of
% 10 log10(power) over 'window' minus its mean over 'baseline' in the same
% epoch (dB); the statistic is the mean of these changes over the band's
% frequencies. Each frequency counts equally, whatever its power (power falls
% steeply with frequency), and a frequency missing at some times (windows that
% grow with 'cycles', blanked samples) does not shift the average. The log
% makes the epoch differences close to symmetric. The test is a one-sample t
% over epochs, with a sign-flip permutation null (the epoch is the unit: each
% epoch's window and baseline are exchanged). For each contact, all 2^K sign
% patterns of its K epochs are used when that is no more than n_perm + 1,
% random ones otherwise (one set for every such contact). On CCEP data each
% contact uses only the epochs it averages in S (its own stimulation left out).
%
% Inputs:
%   S    output of ieeglab_spectral, type 'timefreq', with 'keeptrials', true
%        or 'band'
%
% Options (name/value):
%   'band'      [low high] Hz. Default S.band (70-150 Hz).
%   'window'    [start end] ms tested. Default S.band_window.
%   'baseline'  [start end] ms. Default S.baseline. As in ieeglab_spectral, only
%               times whose whole window lies inside it are used, per frequency.
%   'tail'      'both' (default), 'positive' (increase) or 'negative'.
%   'n_perm'    random sign patterns when not all are used. Default 10000.
%   'correct'   'fdr' (Benjamini-Hochberg, default), 'by', 'bonferroni' or 'none',
%               across contacts (across contacts and times for 'per_time').
%   'alpha'     default 0.05
%   'min_trials' contacts with fewer epochs are not tested. Default 5.
%   'per_time'  also test every time of 'window' against the baseline (as
%               FieldTrip's activation-versus-baseline statistic). Default false.
%               The power at one time is more skewed than its mean over the
%               baseline, so a one-sided test per time can be liberal (on null
%               noise with 5 frequencies: 6.6% below p 0.05 for 'positive', 3.9%
%               for 'negative'; 'both': 5.3% below 0.05, 1.3% below 0.01);
%               prefer 'both' or the window test.
%   'seed'      random sign patterns. Default 1 (same result every run).
%   'verbose'   default true
%
% Output T (one row per contact): channel, band_db (mean over epochs of the
% statistic), t, p, p_adj, significant, n_trials, p_min (the smallest p
% possible: 2/2^K exact two-sided, 1/2^K one-sided, 1/(n_perm + 1) sampled),
% testable (false when it cannot reach significance; see ieeglab_testable).
% Bad contacts and contacts with too few epochs have NaN p and are not counted
% in the correction.
% Output R: the options, times, band_freqs, baseline_mask (band frequencies x
% times), baseline_times (used by at least one frequency), window_times, the
% epochs' statistic d (contacts x epochs, NaN where unused), exact and
% n_patterns per contact, and with 'per_time' t_tc, p_tc, p_adj_tc, sig_tc
% (contacts x times of R.times_tc).
%
% Cedric Cannard, iEEGLAB, 2026

p = inputParser;
p.addParameter('band', []);
p.addParameter('window', []);
p.addParameter('baseline', []);
p.addParameter('tail', 'both', @(x) any(strcmpi(x, {'both', 'positive', 'negative'})));
p.addParameter('n_perm', 10000);
p.addParameter('correct', 'fdr', @(x) any(strcmpi(x, {'fdr', 'by', 'bonferroni', 'none'})));
p.addParameter('alpha', 0.05);
p.addParameter('min_trials', 5);
p.addParameter('per_time', false);
p.addParameter('seed', 1);
p.addParameter('verbose', true);
p.parse(varargin{:});
o = p.Results;
o.tail = lower(char(o.tail)); o.correct = lower(char(o.correct));

if ~isfield(S, 'trials') || ~isfield(S, 'times')
    error('ieeglab_band_stats:needsTrials', ...
        'S must come from ieeglab_spectral(EEG, ''type'', ''timefreq'', ''keeptrials'', true).');
end
if isempty(o.band),     o.band = S.band; end
if isempty(o.window),   o.window = S.band_window; end
if isempty(o.baseline), o.baseline = S.baseline; end

[C, ~, Nt, N] = size(S.trials, 1:4);
toi = S.times;
tf = S.freqs;                                    % frequencies of S.trials ('keeptrials', 'band': the band's)
if isfield(S, 'trials_freqs'), tf = S.trials_freqs; end
inB = tf >= o.band(1) & tf <= o.band(2);
if ~any(inB)
    error('ieeglab_band_stats:bandOutside', 'No frequency of the single-trial power (%s Hz) in the band %s Hz.', ...
        local_range(tf), mat2str(o.band));
end
if any(S.freqs >= o.band(1) & S.freqs <= o.band(2) & ~ismember(S.freqs, tf))
    error('ieeglab_band_stats:bandNotKept', 'S.trials holds %s Hz only, not the whole band %s Hz.', ...
        local_range(tf), mat2str(o.band));
end
fB = tf(inB);                                              % band frequencies of S.trials
[~, loc] = ismember(fB, S.freqs);
ws = S.window_s; if isscalar(ws), ws = repmat(ws, size(S.freqs)); end
if isfield(S, 'srate')
    hMs = (ws(loc) * S.srate - 1) / 2 / S.srate * 1000;     % half-window, ms, as ieeglab_spectral
else
    hMs = ws(loc) * 1000 / 2;
end
inBaseF = toi >= o.baseline(1) + hMs(:) & toi <= o.baseline(2) - hMs(:);   % (band freqs x times)
inWin = toi >= o.window(1) & toi <= o.window(2);
if ~any(inBaseF, 'all')
    error('ieeglab_band_stats:noBaseline', 'No window of the band fits inside the baseline %s ms.', mat2str(o.baseline));
end
noBase = ~any(inBaseF, 2);
if any(noBase) && o.verbose
    warning('ieeglab_band_stats:noBaselineFreq', 'No window fits inside the baseline %s ms at %s Hz: left out.', ...
        mat2str(o.baseline), mat2str(fB(noBase)));
end
if ~any(inWin), error('ieeglab_band_stats:noWindow', 'No time of the window %s ms.', mat2str(o.window)); end

% ---------------- statistic of each epoch ----------------
% Per frequency: window mean minus baseline mean of 10 log10(power), then the
% mean over frequencies (a frequency missing at some times cannot bias it).
nF = numel(fB); tw = find(inWin); nW = numel(tw);
d = nan(C, N);                                              % (contacts x epochs), dB
if o.per_time, Dt = nan(C, nW, N); end                      % (contacts x window times x epochs)
blk = max(1, floor(2^27 / (nF * Nt * N * 8)));              % ~128 MB of double per block of contacts
for c0 = 1:blk:C
    cc = c0:min(C, c0 + blk - 1); nc = numel(cc);
    P = double(S.trials(cc, inB, :, :));                    % (nc x nF x times x epochs)
    P(P <= 0) = NaN;
    L = 10 * log10(P);
    bF = nan(nc, nF, 1, N);
    for k = find(~noBase(:))'
        bF(:, k, 1, :) = mean(L(:, k, inBaseF(k, :), :), 3, 'omitnan');
    end
    dW = L(:, :, tw, :) - bF;                               % (nc x nF x nW x N)
    d(cc, :) = reshape(mean(mean(dW, 3, 'omitnan'), 2, 'omitnan'), nc, N);
    if o.per_time, Dt(cc, :, :) = reshape(mean(dW, 2, 'omitnan'), nc, nW, N); end
end
use = true(C, N);
if isfield(S, 'use'), use = S.use; end
use = use & isfinite(d);
use(S.bad, :) = false;
d(~use) = NaN;
K = sum(use, 2);
tested = K >= o.min_trials;

% ---------------- sign-flip test ----------------
[tObs, pv, exact, J] = local_test(d, use, o.n_perm, o.seed, o.tail);
pv(~tested) = NaN; tObs(~tested) = NaN;

ok = ~isnan(pv);
padj = nan(C, 1);
padj(ok) = local_correct(pv(ok), o.correct);
sig = padj < o.alpha;
tied = 1 + strcmp(o.tail, 'both');                            % s and -s give the same |t|
pMin = nan(C, 1);
pMin(tested & exact) = min(1, tied ./ 2.^K(tested & exact));
pMin(tested & ~exact) = 1 / (o.n_perm + 1);

T = table(S.labels(:), mean(d, 2, 'omitnan'), tObs, pv, padj, sig, K, pMin, true(C, 1), ...
    'VariableNames', {'channel', 'band_db', 't', 'p', 'p_adj', 'significant', 'n_trials', 'p_min', 'testable'});
T.band_db(~tested) = NaN;
if any(ok)
    [testable, msg] = ieeglab_testable(pv(ok), pMin(ok), ones(nnz(ok), 1), @(x) local_correct(x, o.correct), o.alpha);
    T.testable(ok) = testable;
    if ~isempty(msg), warning('ieeglab_band_stats:untestable', '%s', msg); end
end

R = struct('opt', o, 'times', toi, 'band_freqs', fB, 'baseline_mask', inBaseF, ...
    'baseline_times', toi(any(inBaseF, 1)), 'window_times', toi(inWin), 'd', d, 'exact', exact, 'n_patterns', J);

% ---------------- every time of the window against the baseline ----------------
if o.per_time
    ut = reshape(use, C, 1, N) & isfinite(Dt);
    Kt = sum(ut, 3);
    [tt, pt] = local_test(reshape(Dt, C * nW, N), reshape(ut, C * nW, N), o.n_perm, o.seed, o.tail);
    keep = repmat(tested, nW, 1) & Kt(:) >= o.min_trials;
    pt(~keep) = NaN; tt(~keep) = NaN;
    pa = nan(size(pt));
    pa(~isnan(pt)) = local_correct(pt(~isnan(pt)), o.correct);
    R.times_tc = toi(tw);
    R.t_tc = reshape(tt, C, numel(tw));
    R.p_tc = reshape(pt, C, numel(tw));
    R.p_adj_tc = reshape(pa, C, numel(tw));
    R.sig_tc = R.p_adj_tc < o.alpha;
end

if o.verbose
    if ~any(ok), pat = 'no contact tested';
    elseif all(exact(ok)), pat = sprintf('all sign patterns (up to %d)', max(J(ok)));
    elseif ~any(exact(ok)), pat = sprintf('%d random sign patterns', o.n_perm);
    else, pat = sprintf('all or %d random sign patterns', o.n_perm);
    end
    fprintf('[band stats] %g-%g Hz, %s ms vs %s ms: %d of %d contacts differ from their baseline (%s, %s, %s, alpha %g).\n', ...
        o.band(1), o.band(2), mat2str(o.window), mat2str(o.baseline), nnz(sig), nnz(ok), o.tail, pat, o.correct, o.alpha);
end
end

% ===================== local helpers =====================

function [tObs, p, exact, J] = local_test(D, use, nPerm, seed, tail)
% D, use: (cells x epochs). A cell with K used epochs gets all 2^K sign patterns
% of those epochs when 2^K <= nPerm + 1 (cells grouped by K, their used epochs
% moved to the front), otherwise the random patterns over all epochs with the
% unused ones at 0 (one set, the same for every such cell).
[nC, N] = size(D);
D(~use) = 0;
K = sum(use, 2);
tObs = nan(nC, 1); p = nan(nC, 1); J = nan(nC, 1);
exact = 2.^K <= nPerm + 1;
for k = unique(K(exact & K >= 2))'
    r = find(exact & K == k);
    [~, ord] = sort(~use(r, :), 2);                           % used epochs first (stable)
    Dk = reshape(D(sub2ind([nC N], repmat(r, 1, k), ord(:, 1:k))), numel(r), k);
    Sg = double(dec2bin(0:2^k-1, k) == '0')' * 2 - 1;          % (k x 2^k), first column all +1
    [tObs(r), p(r)] = local_perm(Dk, repmat(k, numel(r), 1), Sg, tail);
    J(r) = 2^k;
end
rnd = ~exact;
if any(rnd)
    rs = RandStream('mt19937ar', 'Seed', seed);
    Sg = [ones(N, 1), randi(rs, 2, N, nPerm) * 2 - 3];
    [tObs(rnd), p(rnd)] = local_perm(D(rnd, :), K(rnd), Sg, tail);
    J(rnd) = nPerm + 1;
end
end

function [tObs, p] = local_perm(D, K, Sg, tail)
% D: (cells x epochs), 0 where unused; K: (cells x 1) epochs used; Sg: (epochs x
% patterns), column 1 the observed. t for every pattern at once: sign flips leave
% sum(d.^2) unchanged, so each pattern's variance is (Q - K m^2) / (K - 1) from its
% mean m (exact, in double). Patterns in blocks, cells x patterns near 64 MB.
nC = size(D, 1);
Q = sum(D.^2, 2);
J = size(Sg, 2);
blk = max(1, floor(2^23 / max(1, nC)));
tObs = local_t(D * Sg(:, 1), Q, K);
obs = local_side(tObs, tail);
thr = obs - 1e-12 * abs(obs);                                 % ties within rounding count
thr(obs == Inf) = Inf;                                        % (Inf - Inf would be NaN)
cnt = zeros(nC, 1);
for j0 = 1:blk:J
    jj = j0:min(J, j0 + blk - 1);
    tj = local_t(D * Sg(:, jj), Q, K);
    cnt = cnt + sum(local_side(tj, tail) >= thr, 2);
end
% Column 1 is the observed pattern: all 2^K patterns when exact, otherwise
% (1 + exceedances) / (n_perm + 1)
p = cnt / J;
end

function t = local_t(Ssum, Q, K)
m = Ssum ./ K;
v = (Q - K .* m.^2) ./ (K - 1);
v(v < 0) = 0;                                    % rounding when every |d| is equal
t = m ./ sqrt(v ./ K);
t(m == 0) = 0;
end

function s = local_side(t, tail)
switch tail
    case 'positive', s = t;
    case 'negative', s = -t;
    otherwise,       s = abs(t);
end
s(isnan(s)) = -Inf;
end

function padj = local_correct(p, method)
n = numel(p);
padj = p;
switch method
    case 'bonferroni'
        padj = min(1, p * n);
    case {'fdr', 'by'}
        [ps, ord] = sort(p(:));
        q = ps .* n ./ (1:n)';
        if strcmp(method, 'by'), q = q * sum(1 ./ (1:n)); end
        q = min(1, flipud(cummin(flipud(q))));
        padj(ord) = q;
end
end

function s = local_ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = local_range(f)
if isempty(f), s = 'none'; else, s = sprintf('%g-%g', min(f), max(f)); end
end
